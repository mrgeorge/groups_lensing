pro doall_xgroups_matt_paper,minRadiusKpc,maxRadiusKpc,nRadiusBins,emp_var=emp_var

; outline:
; measure lensing signal around full stacks with different centers
; measure lensing signal around subset stacks with discrepant centers
; model signal from full stack with best center
; model signal from subset stacks with best centers and compare (sensitive to removing high M, low z, etc?)
; construct offset model for full stacks and compare with measured signal
; construct offset model for subset stacks with wrong centers vs. full or subset model of good center

box_factor = 20

; Set paths for input files
path='/users/alexie/Work/Weak_lensing/'
infile_source  = path+'GG_cat_2006/gglensing_source_v1.7.fits'             ; Using the new catalog (photoz version 1.7)
;infile_lens    = path+'Group_cat_june_2010/group7.fits'
;infile_source = '/Users/alexie/Work/GroupCatalogs/lensing15_20101005.fits'
;infile_lens = '/Users/alexie/Work/GroupCatalogs/cosmos_xgroups_20101005.fits'
;infile_lens = '/Users/mgeorge/data/cosmos/code/group5_20110114.fits'
;infile_lens = '/Users/mgeorge/data/cosmos/code/group5_20110201match.fits'
;infile_lens = '/Users/mgeorge/data/cosmos/code/group5_20101005match.fits'
;infile_lens = '/Users/mgeorge/data/cosmos/code/group5_20110114matchbox.fits'
infile_lens = '/Users/mgeorge/data/cosmos/code/group6_20110114.fits'

; Set paths for output files
dirName='bin_'+string(minRadiusKpc,format='(I0)')+'_'+string(maxRadiusKpc,format='(I0)')+'_'+string(nRadiusBins,format='(I0)')
if(keyword_set(emp_var)) then dirName += '_emp'
home=getenv('HOME')
mgPath=home+'/data/cosmos/lensing/test0201f/'+dirName+'/'
if(NOT(file_test(mgPath))) then begin
   file_mkdir,mgPath
   file_mkdir,mgPath+'/plots'
endif else begin
   print, 'Overwrite to '+mgPath+' ?'
;   stop
endelse
dir            = mgPath

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

allnames=['mmgg_scale','mmgg_r200','mmgg2_r200','xray','cm','cm_iter','cl','mlgg_scale','mlgg_r200']
ptSrcAll=[2,2,2,0,0,0,0,2,2]
cen=['mmgg_r200','mmgg2_r200','xray','cm','cm_iter','cl','mlgg_r200']
ptSrcCen=[2,2,0,0,0,0,2]
ref=replicate('mmgg_scale',n_elements(cen))
ptSrcRef=replicate(2,n_elements(cen))

; For one center:
for i=0, n_elements(allnames)-1 do begin
    print,'---------------------------'
    print,allnames[i]
    lens_file = strcompress(dir+'center_'+allnames[i]+'.fits',/remove_all)
    plot_file = strcompress(dir+'Plots/center_'+allnames[i],/remove_all)
    fit_t[0]=ptSrcAll[i]  ; set to include or exclude point source from fit
;    run_gg_offset, infile_source, infile_lens, file, 40, 12344, [35, 50, 0.0, 1.5], box_factor,2, /xgroups,/usespecz,/stackx,center=allnames[i]
;    plot_halofit, [0.5,0.9],file, plot,fit_t, this_m200, this_err_m200, /makeplot,/quick_c,/groups,/use_200,/stackx,center=allnames[i]
    run_gg_offset, infile_source, infile_lens, lens_file, 40, minRadiusKpc, maxRadiusKpc, nRadiusBins, [35, 50, 0.0, 1.5], box_factor,2, /xgroups,/usespecz,center=allnames[i],/stackx,emp_var=keyword_set(emp_var)

   ; Fit model to the lensing signal
   run_ds_mcmc, lens_file, fit_t, rob_p_mean, rob_p_sigma
   ; Plot the results, with and without the models
   plot_lensing_results,lens_file,plot_file,rob_p_mean,fit_t,/stackx,/use_m200
   plot_lensing_results,lens_file,plot_file+'_models',rob_p_mean,fit_t,/stackx,/use_m200,/models
endfor

; Comparing centers:
for i=0,n_elements(cen)-1 do begin

   ; Start by measuring and modeling the signal around the good center
   print,'---------------------------'
   print,ref[i],' and  ',cen[i]
   lens_file = strcompress(dir+'center_'+ref[i]+'_'+cen[i]+'.fits',/remove_all)
   plot_file = strcompress(dir+'plots/center_'+ref[i]+'_'+cen[i],/remove_all)
   run_gg_offset, infile_source, infile_lens, lens_file, 40, minRadiusKpc, maxRadiusKpc, nRadiusBins, [35, 50, 0.0, 1.5], box_factor,2, /xgroups,/usespecz,center=ref[i],refcen=cen[i],/stackx,emp_var=keyword_set(emp_var)

   ; Fit model to the lensing signal around the good center
   fit_t[0]=ptSrcRef[i]            ; set to include or exclude point source from model
   run_ds_mcmc, lens_file, fit_t, rob_p_mean, rob_p_sigma
   ; Plot the results, with and without the models
   plot_lensing_results,lens_file,plot_file,rob_p_mean,fit_t,/stackx,/use_m200
   plot_lensing_results,lens_file,plot_file+'_models',rob_p_mean,fit_t,/stackx,/use_m200,/models


  ; Now measure the lensing signal around the other center and plot both models
   print,'---------------------------'
   print,cen[i],' and  ',ref[i]
   lens_file = strcompress(dir+'center_'+cen[i]+'_'+ref[i]+'.fits',/remove_all)
   plot_file = strcompress(dir+'plots/center_'+cen[i]+'_'+ref[i],/remove_all)
   run_gg_offset, infile_source, infile_lens, lens_file, 40, minRadiusKpc, maxRadiusKpc, nRadiusBins, [35, 50, 0.0, 1.5], box_factor,2, /xgroups,/usespecz,center=cen[i],refcen=ref[i],/stackx,emp_var=keyword_set(emp_var)

   ; Use the results from the good center model fit (rob_p_mean)
   ;  rather than rerunning run_ds_mcmc
   fit_t[0]=ptSrcCen[i]            ; set to include or exclude point source from model
   ; Plot the results with and without models
   plot_lensing_results,lens_file,plot_file,rob_p_mean,fit_t,/stackx,/use_m200,center=cen[i],refcen=ref[i]
   plot_lensing_results,lens_file,plot_file+'_models',rob_p_mean,fit_t,/stackx,/use_m200,center=cen[i],refcen=ref[i],/models

endfor


;-------------------------------
; CONCENTRATION
;-------------------------------

fit_t2 = [$
2,$             ; 0  M0    : baryonic mass
1,$             ; 1  M_vir : NFW Mass
1,$             ; 2  C     : NFW concentration
0,$             ; 5  alpha : fraction
0,$             ; 6  bias
0 ]             ; 6  m_sigma

m200        = fltarr(3)
err_m200    = fltarr(3)
c200        = fltarr(3)
err_c200    = fltarr(3)

center_type        = 4
second_center_type = -1

dir            = path+'Results_group/Conc/'

;file = dir+'alexie_bcg_v1.fits'
;plot = dir+'Plots/alexie_bcg_v1'
;run_gg, infile_source, infile_lens, file, 40, 1235, [30, 42.3, 0.2, 0.9], box_factor,2, /xgroups,/usespecz,typecenter=100    
;plot_halofit, [0.5,0.9],file, plot,fit_t2, this_m200, this_err_m200, p_mean, p_sigma, /makeplot,/quick_c,/groups,/use_200
;m200[0]        = p_mean[0]
;err_m200[0]    = p_sigma[0]
;c200[0]        = p_mean[1]
;err_c200[0]    = p_sigma[1]

file = dir+'alexie_bcg_v2.fits'
plot = dir+'Plots/alexie_bcg_v2'
;run_gg, infile_source, infile_lens, file, 40, 1235, [30, 42.5, 0.2, 0.9], box_factor,2, /xgroups,/usespecz,typecenter=100,/stackx
;just_plot, file, plot, stackx=1
;plot_halofit, [0.5,0.9],file, plot,fit_t2, this_m200, this_err_m200, p_mean, p_sigma, /makeplot,/quick_c,/groups,/use_200
;m200[1]        = p_mean[0]
;err_m200[1]    = p_sigma[0]
;c200[1]        = p_mean[1]
;err_c200[1]    = p_sigma[1]

file = dir+'alexie_bcg_v3.fits'
plot = dir+'Plots/alexie_bcg_v3'
;run_gg, infile_source, infile_lens, file, 40, 1235, [42.5, 50, 0.2, 0.9], box_factor,2, /xgroups,/usespecz,typecenter=100,/stackx 
;just_plot, file, plot, stackx=1
;plot_halofit, [0.5,0.9],file, plot,fit_t2, this_m200, this_err_m200, p_mean, p_sigma, /makeplot,/quick_c,/groups,/use_200
;m200[2]        = p_mean[0]
;err_m200[2]    = p_sigma[0]
;c200[2]        = p_mean[1]
;err_c200[2]    = p_sigma[1]

; NOTE: Conc here is in log10 !!
;save, m200, filename='/Users/alexie/Work/Weak_lensing/Results_group/Conc/m200_v1_conc.sav'
;save, err_m200, filename='/Users/alexie/Work/Weak_lensing/Results_group/Conc/err_m200_v1_conc.sav'
;save, c200, filename='/Users/alexie/Work/Weak_lensing/Results_group/Conc/c200_v1_conc.sav'
;save, err_c200, filename='/Users/alexie/Work/Weak_lensing/Results_group/Conc/err_c200_v1_conc.sav'

;compare_conc
end
