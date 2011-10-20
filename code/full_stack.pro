pro model_full_stack,cenNames,ptSrc,infile_source,infile_lens,lensOutFileArr,$
                     innerRadiusKpc,secondRadiusKpc,maxRadiusKpc,nRadiusBins,$
                     minLensZ,maxLensZ,minLensMass,maxLensMass,$
                     stackx=stackx,emp_var=emp_var

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
0,$             ; 3  alpha : fraction
0,$             ; 4  bias
0,$             ; 5  m_sigma
0]              ; 6  offset
;----------------------------------------------------------------------
; Various tests 
; zscheme 0 -> 68 % cut + sigma cut=0.1
; zscheme 1 -> 68 % cut + sigma cut=0.2
; zscheme 2 -> remove double peak objects + 68 % cut + sigma cut=0.1
; zscheme 3 -> remove double peak objects + 68 % cut + sigma cut=0.2
;----------------------------------------------------------------------
zscheme=2 ; from ALs code
box_factor=20 ; from ALs code

; create 2d array to save fit types, 1 row for each center
fitTypeAll=rebin(fit_t,n_elements(fit_t),n_elements(cenNames))
fitTypeAll[0,*]=ptSrc

; repeat for offset model
fitTypeAll2=fitTypeAll
fitTypeAll2[6,*]=1 ; this leaves the offset as a free parameter

; For one center:
for i=0, n_elements(cenNames)-1 do begin
    print,'---------------------------'
    print,cenNames[i]

    run_gg_offset, infile_source, infile_lens, lensOutFileArr[i], innerRadiusKpc, secondRadiusKpc, maxRadiusKpc, nRadiusBins, minLensZ, maxLensZ, minLensMass, maxLensMass, box_factor, zscheme, /xgroups,/usespecz,center=cenNames[i],stackx=keyword_set(stackx),emp_var=keyword_set(emp_var)

   ; Fit model to the lensing signal
    ; fit pars are returned as rob_p_mean and also written to lensOutFile
   run_ds_mcmc, lensOutFileArr[i], fitTypeAll[*,i], rob_p_mean, rob_p_sigma, stackx=keyword_set(stackx)

   ; Fit offset model to the lensing signal
    ; fit pars are returned as rob_p_mean and also written to lensOutFile
   run_ds_mcmc, lensOutFileArr[i], fitTypeAll2[*,i], rob_p_mean2, rob_p_sigma2, stackx=keyword_set(stackx),/fast

   ; Plot the results, with and without the models
;   plot_lensing_results,lensOutFileArr[i],plotFileArr[i],rob_p_mean,fitTypeAll[*,i],stackx=keyword_set(stackx),/use_m200
;   plot_lensing_results,lensOutFileArr[i],plotFileArr[i]+'_models',rob_p_mean,fitTypeAll[*,i],stackx=keyword_set(stackx),/use_m200,/models
endfor

end

pro table_full_stack,cenTitles,lensFileArr,tableFile,stackx=stackx
  
; reorder center arrays so that table order corresponds with plot from
; top to bottom
reorder=[6,7,4,5,2,3,0,1]
cenTitles=cenTitles[reorder]
lensFileArr=lensFileArr[reorder]

openw,u,tableFile,/get_lun,width=1000
printf,u,'# Center  Mcen  Mnfw  dMnfw  chisq  Mnfw_off  dMnfw_off  Roff  dRoff  chisq_off'

for ii=0,n_elements(cenTitles)-1 do begin
   str=mrdfits(lensFileArr[ii],1)

   ; restrict to points with enough sources
   sel=where(str.e1_num GE 10 AND str.plot_radius_kpc GT 10)
   if(keyword_set(stackx)) then begin
      stop ; need to rewrite this to determine rnfw if needed
      x    = str.plot_radius_kpc[sel]*rnfw
      y    = str.we1_mean[sel]
      yerr = str.we1_error[sel]
   endif else begin
      x    = str.plot_radius_kpc[sel]/1e3
      y    = str.we1_mean[sel]
      yerr = str.we1_error[sel]
   endelse
   
   chisq=get_ds_chisq(str.fit_type,str.p_mean,str,x,y,yerr,dof=dof,/use_m200)
   chisq2=get_ds_chisq(str.fit_type2,str.p_mean2,str,x,y,yerr,dof=dof2,/use_m200)

   if(str.msun_lens EQ 0.) then mcen='-' else mcen=string(str.msun_lens,format='(F4.1)')
      
   if(ii EQ 0) then printf,u,'# dof (centered, offset): ',dof,dof2

   printf,u,cenTitles[ii],mcen,str.p_mean,str.p_sigma,chisq,str.p_mean2[0],str.p_sigma2[0],1000.*str.p_mean2[1],1000.*str.p_sigma2[1],chisq2, $
          FORMAT='(A15," &",A5," &",F6.2," $\pm$",F6.2," &",F5.1," &",F6.2," $\pm$",F6.2," &",F6.1," $\pm$",F6.1," &",F5.1," \\")'
endfor

close,u
free_lun,u
end

pro full_stack,innerRadiusKpc,secondRadiusKpc,maxRadiusKpc,nRadiusBins,stackx=stackx,emp_var=emp_var,loz=loz,hiz=hiz,loM=loM,hiM=hiM

; measure lensing signal for the full stack of groups around different centers
; plot results in one big figure for paper

; code evolved from doall_xgroups_matt_paper.pro
; split into this script for single full stacks and a separate script
; diff_stack.pro to measuring lensing signal comparing two different
; centers

; sample call: full_stack,10,50,2000,8,/stackx,/emp_var

if(keyword_set(loz)) then begin
   minLensZ=0.0
   maxLensZ=0.5
   zExt='_loz'
endif else if(keyword_set(hiz)) then begin
   minLensZ=0.5
   maxLensZ=1.0
   zExt='_hiz'
endif else begin
   minLensZ=0.0
   maxLensZ=1.0
   zExt=''
endelse
if(keyword_set(loM)) then begin
   minLensMass=12.0
   maxLensMass=13.5
   mExt='_loM'
endif else if(keyword_set(hiM)) then begin
   minLensMass=13.5
   maxLensMass=15.0
   mExt='_hiM'
endif else begin
   minLensMass=12.0
   maxLensMass=15.0
   mExt=''
endelse
if(keyword_set(stackx)) then sxExt='_sx' else sxExt=''
if(keyword_set(emp_var)) then empExt='_emp' else empExt=''

; Set paths for input files
infile_source='/Users/alexie/Work/Weak_lensing/GG_cat_2006/gglensing_source_v1.7.fits' ; Using the new catalog (photoz version 1.7)
infile_lens = '/Users/alexie/Work/GroupCatalogs/cosmos_xgroups_20110209.fits' ; group catalog with centers

; Set paths for output files
dirName='bin_'+string(innerRadiusKpc,format='(I0)')+'_'+string(secondRadiusKpc,format='(I0)')+'_'+string(maxRadiusKpc,format='(I0)')+'_'+string(nRadiusBins,format='(I0)')+sxExt+empExt+zExt+mExt
fileDir='~/data/cosmos/groups_lensing/outfiles/'+dirName+'/'
plotDir='~/data/cosmos/groups_lensing/plots/'+dirName+'/'
if(NOT(file_test(fileDir))) then file_mkdir,fileDir
if(NOT(file_test(plotDir))) then file_mkdir,plotDir

; Centers ordered bottom left to top right on plot
cenNames=['cl','xray','cn','cm','mlgg_scale','mlgg_r200','mmgg_scale','mmgg_r200']
cenTitlesTex=['CL','X-ray','CN','CM','MLGG_{scale}','MLGG_{R200}','MMGG_{scale}','MMGG_{R200}']
cenTitles=textoidl(cenTitlesTex)
ptSrc=[0,0,0,0,2,2,2,2] ; for fit_t

lensOutFileArr=strcompress(fileDir+'center_'+cenNames+'.fits',/remove_all)
plotFileArr=strcompress(plotDir+'center_'+cenNames,/remove_all)
fullPlotFile=plotDir+'full_stacks.eps'
tableFile=plotDir+'full_stack_fit_pars.tex'

; Measure and Model the lensing signal
;model_full_stack,cenNames,ptSrc,infile_source,infile_lens,lensOutFileArr,$
;                     innerRadiusKpc,secondRadiusKpc,maxRadiusKpc,nRadiusBins,$
;                     minLensZ,maxLensZ,minLensMass,maxLensMass,$
;                     stackx=stackx,emp_var=emp_var

; Make multi-panel plot to compare full stacks
plot_full_stacks,cenTitles,lensOutFileArr,fullPlotFile,stackx=keyword_set(stackx),/use_m200

; Print latex formatted table with fit parameters
table_full_stack,cenTitlesTex,lensOutFileArr,tableFile,stackx=stackx

end
