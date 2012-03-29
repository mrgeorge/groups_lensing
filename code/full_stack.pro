pro model_full_stack,cenNames,ptSrc,infile_source,infile_lens,lensOutFileArr,$
                     innerRadiusKpc,secondRadiusKpc,maxRadiusKpc,nRadiusBins,$
                     minLensZ,maxLensZ,minLensMass,maxLensMass,$
                     stackx=stackx,emp_var=emp_var,conc=conc,cenFree=cenFree,subhalo=subhalo

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

if(n_elements(conc) GT 0) then fit_t[2]=1
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
if(n_elements(cenFree) GT 0) then fitTypeAll[0,*]=1 $
else fitTypeAll[0,*]=ptSrc


; repeat for offset model
fitTypeAll2=fitTypeAll
fitTypeAll2[6,*]=1 ; this leaves the offset as a free parameter

; For one center:
for i=0, n_elements(cenNames)-1 do begin
    print,'---------------------------'
    print,cenNames[i]

    run_gg_offset, infile_source, infile_lens, lensOutFileArr[i], innerRadiusKpc, secondRadiusKpc, maxRadiusKpc, nRadiusBins, minLensZ, maxLensZ, minLensMass, maxLensMass, box_factor, zscheme, /xgroups,/usespecz,center=cenNames[i],stackx=keyword_set(stackx),emp_var=keyword_set(emp_var),subhalo=subhalo

   ; Fit centered model to the lensing signal
    ; fit pars are returned as rob_p_mean and also written to lensOutFile
   run_ds_mcmc, lensOutFileArr[i], fitTypeAll[*,i], rob_p_mean, rob_p_sigma, stackx=keyword_set(stackx),/fast,/ps

   ; Fit offset model to the lensing signal
    ; fit pars are returned as rob_p_mean and also written to lensOutFile
   run_ds_mcmc, lensOutFileArr[i], fitTypeAll2[*,i], rob_p_mean2, rob_p_sigma2, stackx=keyword_set(stackx),/fast,/ps,/off3dMax

   ; Plot the results, with and without the models
;   plot_lensing_results,lensOutFileArr[i],plotFileArr[i],rob_p_mean,fitTypeAll[*,i],stackx=keyword_set(stackx),/use_m200
;   plot_lensing_results,lensOutFileArr[i],plotFileArr[i]+'_models',rob_p_mean,fitTypeAll[*,i],stackx=keyword_set(stackx),/use_m200,/models
endfor

end

pro table_full_stack,cenTitles,lensFileArr,tableFile,stackx=stackx,conc=conc,cenFree=cenFree
  
; reorder center arrays so that table order corresponds with plot from
; top to bottom
reorder=[6,7,4,5,2,3,0,1]
cenTitles=cenTitles[reorder]
lensFileArr=lensFileArr[reorder]

openw,u,tableFile,/get_lun,width=1000
if(n_elements(conc) EQ 0 AND n_elements(cenFree) EQ 0) then printf,u,'# Center  Mcen  Mnfw  dMnfw  Conc  chisq  Mnfw_off  dMnfw_off  Conc_off  Roff  dRoff  chisq_off' $
else if(n_elements(conc) EQ 0 AND n_elements(cenFree) GT 0) then printf,u,'# Center  Mcen  dMCen  Mnfw  dMnfw  Conc  chisq  Mcen_off  dMcen_off  Mnfw_off  dMnfw_off  Conc_off  Roff  dRoff  chisq_off' $
else if(n_elements(conc) GT 0 AND n_elements(cenFree EQ 0)) then printf,u,'# Center  Mcen  Mnfw  dMnfw  Conc  dConc  chisq  Mnfw_off  dMnfw_off  Conc_off  dConc_ff  Roff  dRoff  chisq_off' $
else printf,u,'# Center  Mcen  dMcen  Mnfw  dMnfw  Conc  dConc  chisq  Mcen_off  dMcen_off  Mnfw_off  dMnfw_off  Conc_off  dConc_ff  Roff  dRoff  chisq_off'

for ii=0,n_elements(cenTitles)-1 do begin
   str=mrdfits(lensFileArr[ii],1)

   ; restrict to points with enough sources
   sel=where(str.e1_num GE 10 AND str.plot_radius_kpc GT 20)
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

   fit_type=str.fit_type
   fit_type2=str.fit_type2
   if(fit_type[0] GT 0) then cen_type=str.cen_type
   if(fit_type[6] GT 0) then off_type=str.off_type
   if(fit_type2[0] GT 0) then cen_type2=str.cen_type2
   if(fit_type2[6] GT 0) then off_type2=str.off_type2
   
   ; get Zhao concentrations if conc was not free
   if(n_elements(conc) EQ 0) then begin
      get_ds_model,fit_type,str.p_mean,str.z_lens,str.msun_lens,x,/use_m200,conc=concVal,cen_type=cen_type,off_type=off_type
      get_ds_model,fit_type2,str.p_mean2,str.z_lens,str.msun_lens,x,/use_m200,conc=concVal2,cen_type=cen_type2,off_type=off_type2
   endif

   ; get chi^2 for each model
   chisq=get_ds_chisq(fit_type,str.p_mean,str.z_lens,str.msun_lens,x,y,yerr,dof=dof,/use_m200,cen_type=cen_type,off_type=off_type)
   chisq2=get_ds_chisq(fit_type2,str.p_mean2,str.z_lens,str.msun_lens,x,y,yerr,dof=dof2,/use_m200,cen_type=cen_type2,off_type=off_type2)

   if(n_elements(cenFree) EQ 0) then begin
      if(str.msun_lens EQ 0.) then mcen='-' else mcen=string(str.msun_lens,format='(F4.1)')
   endif else mcen=string(str.p_mean[0],format='(F4.1)')

   if(ii EQ 0) then printf,u,'# dof (centered, offset): ',dof,dof2

   if(n_elements(conc) EQ 0 AND n_elements(cenFree) EQ 0) then printf,u,cenTitles[ii],mcen,str.p_mean,str.p_sigma,concVal,chisq,str.p_mean2[0],str.p_sigma2[0],concVal2,1000.*str.p_mean2[1],1000.*str.p_sigma2[1],chisq2, $
          FORMAT='(A15," &",A5," &",F6.2," $\pm$",F6.2," &",F5.1," &",F5.1," & & ",F6.2," $\pm$",F6.2," &",F5.1," &",F6.1," $\pm$",F6.1," &",F5.1," \\")' $

   else if(n_elements(conc) GT 0 AND n_elements(cenFree) EQ 0) then printf,u,cenTitles[ii],mcen,str.p_mean[0],str.p_sigma[0],str.p_mean[1],str.p_sigma[1],chisq,str.p_mean2[0],str.p_sigma2[0],str.p_mean2[1],str.p_sigma2[1],1000.*str.p_mean2[2],1000.*str.p_sigma2[2],chisq2, $
          FORMAT='(A15," &",A5," &",F6.2," $\pm$",F6.2," &",F5.1," $\pm$",F5.1," &",F5.1," & & ",F6.2," $\pm$",F6.2," &",F5.1," $\pm$",F5.1," &",F6.1," $\pm$",F6.1," &",F5.1," \\")' $

   else if(n_elements(conc) EQ 0 AND n_elements(cenFree) GT 0) then printf,u,cenTitles[ii],mcen,str.p_sigma[0],str.p_mean[1],str.p_sigma[1],concVal,chisq,str.p_mean2[0],str.p_sigma2[0],str.p_mean2[1],str.p_sigma2[1],concVal2,1000.*str.p_mean2[2],1000.*str.p_sigma2[2],chisq2, $
          FORMAT='(A15," &",A5," $\pm$",F4.1," &",F6.2," $\pm$",F6.2," &",F5.1," &",F5.1," & & ",F4.1," $\pm$",F4.1," &",F6.2," $\pm$",F6.2," &",F5.1," &",F6.1," $\pm$",F6.1," &",F5.1," \\")' $

   else if(n_elements(conc) GT 0 AND n_elements(cenFree) GT 0) then printf,u,cenTitles[ii],mcen,str.p_sigma[0],str.p_mean[1],str.p_sigma[1],str.p_mean[2],str.p_sigma[2],chisq,str.p_mean2[0],str.p_sigma2[0],str.p_mean2[1],str.p_sigma2[1],str.p_mean2[2],str.p_sigma2[2],1000.*str.p_mean2[3],1000.*str.p_sigma2[3],chisq2, $
          FORMAT='(A15," &",A5," $\pm$",F4.1," &",F6.2," $\pm$",F6.2," &",F5.1," $\pm$",F5.1," &",F5.1," & & ",F4.1," $\pm$",F4.1," &",F6.2," $\pm$",F6.2," &",F5.1," $\pm$",F5.1," &",F6.1," $\pm$",F6.1," &",F5.1," \\")'
endfor

close,u
free_lun,u
end

pro full_stack,innerRadiusKpc,secondRadiusKpc,maxRadiusKpc,nRadiusBins,stackx=stackx,emp_var=emp_var,loz=loz,hiz=hiz,loM=loM,hiM=hiM,conc=conc,cenFree=cenFree,subhalo=subhalo

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
if(keyword_set(conc)) then concExt='_conc' else concExt=''
if(keyword_set(cenFree)) then cenExt='_cen' else cenExt=''
if(keyword_set(subhalo)) then subExt='_sub'+string(subhalo,format='(F04.1)') else subExt=''

; Set paths for input files
infile_source='/Users/alexie/Work/Weak_lensing/GG_cat_2006/gglensing_source_v1.7.fits' ; Using the new catalog (photoz version 1.7)
;infile_lens = '/Users/alexie/Work/GroupCatalogs/cosmos_xgroups_20110209.fits' ; group catalog with centers
infile_lens = '~/data/cosmos/code/group5_20110914.fits' ; group catalog with centers

; Set paths for output files
dirName='bin_'+string(innerRadiusKpc,format='(I0)')+'_'+string(secondRadiusKpc,format='(I0)')+'_'+string(maxRadiusKpc,format='(I0)')+'_'+string(nRadiusBins,format='(I0)')+sxExt+empExt+concExt+cenExt+subExt+zExt+mExt+"_20110914"
fileDir='~/data/cosmos/groups_lensing/outfiles/'+dirName+'/'
plotDir='~/data/cosmos/groups_lensing/plots/'+dirName+'/'
if(NOT(file_test(fileDir))) then file_mkdir,fileDir
if(NOT(file_test(plotDir))) then file_mkdir,plotDir

; Centers ordered bottom left to top right on plot
cenNames=['cn','cm','cf','xray','mmgg_scale','mmgg_r200','bgg_scale','bgg_r200']
cenTitlesTex=['CN','CM','CF','X-ray','MMGG_{scale}','MMGG_{R200}','BGG_{scale}','BGG_{R200}']
cenTitles=textoidl(cenTitlesTex)
ptSrc=[0,0,0,0,2,2,2,2] ; for fit_t

lensOutFileArr=strcompress(fileDir+'center_'+cenNames+'.fits',/remove_all)
plotFileArr=strcompress(plotDir+'center_'+cenNames,/remove_all)
fullPlotFile=plotDir+'full_stacks_split.eps'
tableFile=plotDir+'full_stack_fit_pars.tex'

; Measure and Model the lensing signal
model_full_stack,cenNames,ptSrc,infile_source,infile_lens,lensOutFileArr,$
                     innerRadiusKpc,secondRadiusKpc,maxRadiusKpc,nRadiusBins,$
                     minLensZ,maxLensZ,minLensMass,maxLensMass,$
                     stackx=stackx,emp_var=emp_var,conc=conc,cenFree=cenFree,subhalo=subhalo

; Make multi-panel plot to compare full stacks
plot_full_stacks,cenTitles,lensOutFileArr,fullPlotFile,stackx=keyword_set(stackx),/use_m200

; Print latex formatted table with fit parameters
table_full_stack,cenTitlesTex,lensOutFileArr,tableFile,stackx=stackx,conc=conc,cenFree=cenFree

end
