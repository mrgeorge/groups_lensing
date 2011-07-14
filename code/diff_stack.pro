pro diff_stack,minRadiusKpc,maxRadiusKpc,nRadiusBins,stackx=stackx,emp_var=emp_var

; measure lensing signal for two stacks of groups where centers disagree
; plot results in one big figure for paper

; code evolved from doall_xgroups_matt_paper.pro
; split into this script for comparing stacks and a separate script
; full_stack.pro for single full stacks

; sample call: diff_stack,50,2000,8,/stackx,/emp_var

; Set paths for input files
infile_source='/Users/alexie/Work/Weak_lensing/GG_cat_2006/gglensing_source_v1.7.fits' ; Using the new catalog (photoz version 1.7)
infile_lens = '/Users/alexie/Work/GroupCatalogs/cosmos_xgroups_20110209.fits' ; group catalog with centers
;infile_lens = '/Users/mgeorge/data/cosmos/code/group6_20110712.fits' ; group catalog with red centers

; Set paths for output files
dirName='bin_'+string(minRadiusKpc,format='(I0)')+'_'+string(maxRadiusKpc,format='(I0)')+'_'+string(nRadiusBins,format='(I0)')
if(keyword_set(stackx)) then dirName += '_sx'
if(keyword_set(emp_var)) then dirName += '_emp'
fileDir='~/data/cosmos/groups_lensing/outfiles/'+dirName+'/'
plotDir='~/data/cosmos/groups_lensing/plots/'+dirName+'/'
if(NOT(file_test(fileDir))) then file_mkdir,fileDir
if(NOT(file_test(plotDir))) then file_mkdir,plotDir

;-------------------------------------------------------------------------
; Fitting parameters
; 0  -> Ignore this componant in the fit 
;       If any parm is the componant is set to 0 then the whole comp will be ignored
; 1  -> Free parameter. Fit this from the data.
; 2  -> Fix this parameter from either model or data (e.g stellar mass)
;-------------------------------------------------------------------------
fit_t = [$
2,$             ; 0  M0    : baryonic mass
1,$             ; 1  R_vir : NFW virial mass
2,$             ; 2  C     : NFW concentration
0,$             ; 5  alpha : fraction
0,$             ; 6  bias
0 ]             ; 6  m_sigma
;----------------------------------------------------------------------
; Various tests 
; zscheme 0 -> 68 % cut + sigma cut=0.1
; zscheme 1 -> 68 % cut + sigma cut=0.2
; zscheme 2 -> remove double peak objects + 68 % cut + sigma cut=0.1
; zscheme 3 -> remove double peak objects + 68 % cut + sigma cut=0.2
;----------------------------------------------------------------------
zscheme=2 ; from ALs code
box_factor=20 ; from ALs code
minLensZ=0.0
maxLensZ=1.0
minLensMass=12.
maxLensMass=15.

; ordered bottom left to top on plot
cenNames=['mlgg_scale','mlgg_r200','mmgg_r200','cl','cm','cn','xray']
cenText=textoidl(['MLGG_{scale}','MLGG_{R200}','MMGG_{R200}','CL','CM','CN','X-ray'])
ptSrcCen=[2,2,2,0,0,0,0] ; for fit_t
refNames=replicate('mmgg_scale',n_elements(cenNames)) ; the "good center" to compare with the ones above
refText=replicate(textoidl('MMGG_{scale}'),n_elements(cenNames))
ptSrcRef=replicate(2,n_elements(cenNames))

; filenames convention - first name is the center used for lensing,
;                        2nd is the one used for comparison (needed
;                        b/c only stacking groups where they disagree)
lensOutFileArrCen=strcompress(fileDir+'center_'+cenNames+'_'+refNames+'.fits',/remove_all)
plotFileArrCen=strcompress(plotDir+'center_'+cenNames+'_'+refNames,/remove_all)
lensOutFileArrRef=strcompress(fileDir+'center_'+refNames+'_'+cenNames+'.fits',/remove_all)
plotFileArrRef=strcompress(plotDir+'center_'+refNames+'_'+cenNames,/remove_all)
diffPlotFile=plotDir+'diff_stacks.eps'

; create 2d array to save fit types, 1 row for each center
fitTypeAllCen=rebin(fit_t,n_elements(fit_t),n_elements(cenNames))
fitTypeAllCen[0,*]=ptSrcCen
fitTypeAllRef=rebin(fit_t,n_elements(fit_t),n_elements(refNames))
fitTypeAllRef[0,*]=ptSrcRef

; arrays to record mean and stddev mass from mcmc
; we fit mass on reference center, and just use that model to get
; offset w/o mcmc
massMeanRef=fltarr(n_elements(refNames))
massErrRef=fltarr(n_elements(refNames))

; Comparing centers:
for i=0,n_elements(cenNames)-1 do begin

   ; Start by measuring and modeling the signal around the good center
   print,'---------------------------'
   print,refNames[i],' and  ',cenNames[i]

   run_gg_offset, infile_source, infile_lens, lensOutFileArrRef[i], minRadiusKpc, maxRadiusKpc, nRadiusBins, minLensZ, maxLensZ, minLensMass, maxLensMass, box_factor, zscheme, /xgroups,/usespecz,center=refNames[i],refcen=cenNames[i],stackx=keyword_set(stackx),emp_var=keyword_set(emp_var)

   ; Fit model to the lensing signal around the good center
   run_ds_mcmc, lensOutFileArrRef[i], fitTypeAllRef[*,i], rob_p_mean, rob_p_sigma
   massMeanRef[i]=rob_p_mean
   massErrRef[i]=rob_p_sigma

   ; Plot the results, with and without the models
   plot_lensing_results,lensOutFileArrRef[i],plotFileArrRef[i],rob_p_mean,fitTypeAllRef[*,i],stackx=keyword_set(stackx),/use_m200
   plot_lensing_results,lensOutFileArrRef[i],plotFileArrRef[i]+'_models',rob_p_mean,fitTypeAllRef[*,i],stackx=keyword_set(stackx),/use_m200,/models


  ; Now measure the lensing signal around the other center and plot both models
   print,'---------------------------'
   print,cenNames[i],' and  ',refNames[i]

   run_gg_offset, infile_source, infile_lens, lensOutFileArrCen[i], minRadiusKpc, maxRadiusKpc, nRadiusBins, minLensZ, maxLensZ, minLensMass, maxLensMass, box_factor, zscheme, /xgroups,/usespecz,center=cenNames[i],refcen=refNames[i],stackx=keyword_set(stackx),emp_var=keyword_set(emp_var)

   ; Use the results from the good center model fit (rob_p_mean)
   ;  rather than rerunning run_ds_mcmc

   ; Plot the results with and without models
   plot_lensing_results,lensOutFileArrCen[i],plotFileArrCen[i],rob_p_mean,fitTypeAllCen[*,i],stackx=keyword_set(stackx),/use_m200,center=cenNames[i],refcen=refNames[i],groupFile=infile_lens
   plot_lensing_results,lensOutFileArrCen[i],plotFileArrCen[i]+'_models',rob_p_mean,fitTypeAllCen[*,i],stackx=keyword_set(stackx),/use_m200,center=cenNames[i],refcen=refNames[i],groupFile=infile_lens,/models

endfor

; plot halo masses and uncertainties from different centers for comparison
set_plot,'ps'
simpctable
!p.font=0
!p.thick=3
!x.thick=3
!y.thick=3
!p.charthick=1.2
!p.charsize=1.2
device,filename=plotDir+'mass_comparison_diff.eps',/encapsul,/color,/helvetica,xsize=10,ysize=4,/inches

plot,massMeanRef,psym=4,xticks=n_elements(cenNames)-1,xtickn=cenNames,xminor=1,xstyle=2,yrange=[min(massMeanRef-massErrRef),max(massMeanRef+massErrRef)]
oploterror,massMeanRef,massErrRef,psym=4,errthick=3

device,/close

; record the masses and errors to a text file
openw,u,plotDir+'mass_comparison_diff.txt',/get_lun
printf,u,'# Masses for stacks centered on Center1 when they disagree with Center2'
printf,u,'# Masses and errors from ds_mcmc. M200c, h=0.72'
printf,u,'# Center1  Center2  Mass  MassErr'
for ii=0,n_elements(cenNames)-1 do printf,u,refNames[ii],cenNames[ii],massMeanRef[ii],massErrRef[ii]
close,u
free_lun,u

plot_diff_stacks,cenNames,refNames,cenText,refText,lensOutFileArrCen,lensOutFileArrRef,diffPlotFile,infile_lens,massMeanRef,fitTypeAllCen,fitTypeAllRef,stackx=keyword_set(stackx),/use_m200

end
