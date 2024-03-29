pro jackknife_stack, refName, cenName, hiZ=hiZ

; stack random subsamples of groups to determine the distribution of
; masses obtained. This is a test to see whether groups with
; MMGG_scale != MMGG_r200 have anomalously low masses.

;group=mrdfits("~alexie/Work/GroupCatalogs/cosmos_xgroups_20110209.fits",1)
group=mrdfits("~/data/cosmos/code/group5_20110914.fits",1)
sel=where(group.flag_include EQ 1)
group=group[sel]

if(keyword_set(hiZ)) then begin
   sel=where(group.zphot GT 0.7)
   group=group[sel]
endif

minSep=50. ; kpc
if(refName EQ 'mmgg_scale' AND (cenName EQ 'mlgg_scale' OR cenName EQ 'bgg_scale')) then begin
   bad=where(distance(group.alpha_mmgg_scale,group.delta_mmgg_scale,group.alpha_mlgg_scale,group.delta_mlgg_scale)*3600.*group.lensing_r200_mpc*1000./group.lensing_r200_as GT minSep)
endif else if(refName EQ 'mmgg_scale' AND cenName EQ 'mmgg_r200') then begin
   bad=where(distance(group.alpha_mmgg_scale,group.delta_mmgg_scale,group.alpha_mmgg_r200,group.delta_mmgg_r200)*3600.*group.lensing_r200_mpc*1000./group.lensing_r200_as GT minSep)
endif else if(refName EQ 'mmgg_scale' AND (cenName EQ 'mlgg_r200' OR cenName EQ 'bgg_r200')) then begin
   bad=where(distance(group.alpha_mmgg_scale,group.delta_mmgg_scale,group.alpha_mlgg_r200,group.delta_mlgg_r200)*3600.*group.lensing_r200_mpc*1000./group.lensing_r200_as GT minSep)
endif

LXBad=median(group[bad].lx_scale)
mLXBad=median(group[bad].lensing_m200)
zBad=median(group[bad].zphot)

sampleSize=n_elements(bad)

nSamples=1000

LXSamp=fltarr(nSamples)
mLXSamp=fltarr(nSamples) ; median mass inferred from LX-M relation
massSamp=fltarr(nSamples) ; mean mass from WL
massErrSamp=fltarr(nSamples)
zSamp=fltarr(nSamples)
nDiffSamp=intarr(nSamples) ; number of groups where centers disagree

innerRadiusKpc=20.
secondRadiusKpc=70.
maxRadiusKpc=1000.
nRadiusBins=7
stackx=0
emp_var=1
dirName='bin_'+string(innerRadiusKpc,format='(I0)')+'_'+string(secondRadiusKpc,format='(I0)')+'_'+string(maxRadiusKpc,format='(I0)')+'_'+string(nRadiusBins,format='(I0)')
if(keyword_set(stackx)) then dirName += '_sx'
if(keyword_set(emp_var)) then dirName += '_emp'
fileDir='~/data/cosmos/groups_lensing/outfiles/'+dirName+'_20110914/'
if(keyword_set(hiZ)) then zExt='_hiz' else zExt=''
lensOutFile=strcompress(fileDir+'center_'+refName+'_'+cenName+'_jackknife'+zExt+'.fits',/remove_all)
sampGroupFile=fileDir+'group_'+refName+'_'+cenName+'_jackknife'+zExt+'.fits'
infile_source='/Users/alexie/Work/Weak_lensing/GG_cat_2006/gglensing_source_v1.7.fits' ; Using the new catalog (photoz version 1.7)

zscheme=2 ; from ALs code
box_factor=20 ; from ALs code
minLensZ=0.0
maxLensZ=1.0
minLensMass=12.
maxLensMass=15.
fitType = [$
2,$             ; 0  M0    : baryonic mass
1,$             ; 1  R_vir : NFW virial mass
2,$             ; 2  C     : NFW concentration
0,$             ; 3  alpha : fraction
0,$             ; 4  bias
0,$             ; 5  m_sigma
0 ]             ; 6  offset

openw,u,fileDir+'jackknife_results_20110914_'+refName+'_'+cenName+zExt+'.txt',/get_lun,width=1000
printf,u,'# Jackknife sample size',sampleSize
printf,u,'# Run  nDiff  z  LX  mLX  mWL  sig_mWL'
for ii=0,nSamples-1 do begin
   print,'==========================='
   print,'Sample ',ii,' of ',nSamples
   print,'==========================='

   if(ii EQ 0) then sel=bad $ ; list the "bad" sample first for comparison
   else begin
      ; choose N random indices from Ngroups without replacement 
      rand=randomu(seed,n_elements(group))
      sel=(sort(rand))[0:sampleSize-1]
   endelse

   mwrfits,group[sel],sampGroupFile,/create 

   LXSamp[ii]=median(group[sel].lx_scale)
   mLXSamp[ii]=median(group[sel].lensing_m200)
   zSamp[ii]=median(group[sel].zphot)
   
   if(refName EQ 'mmgg_scale' AND (cenName EQ 'mlgg_scale' OR cenName EQ 'bgg_scale')) then begin
      nDiffSamp[ii]=n_elements(where(group[sel].id_mmgg_scale NE group[sel].id_mlgg_scale))
   endif else if(refName EQ 'mmgg_scale' AND cenName EQ 'mmgg_r200') then begin
      nDiffSamp[ii]=n_elements(where(group[sel].id_mmgg_scale NE group[sel].id_mmgg_r200))
   endif else if(refName EQ 'mmgg_scale' AND (cenName EQ 'mlgg_r200' OR cenName EQ 'bgg_r200')) then begin
      nDiffSamp[ii]=n_elements(where(group[sel].id_mmgg_scale NE group[sel].id_mlgg_r200))
   endif

   run_gg_offset, infile_source, sampGroupFile, lensOutFile, innerRadiusKpc, secondRadiusKpc, maxRadiusKpc, nRadiusBins, minLensZ, maxLensZ, minLensMass, maxLensMass, box_factor, zscheme, /xgroups,/usespecz,center=refName,stackx=keyword_set(stackx),emp_var=keyword_set(emp_var)

   ; Fit model to the lensing signal
   run_ds_mcmc, lensOutFile, fitType, rob_p_mean, rob_p_sigma, /fast, stackx=keyword_set(stackx), /noSave, /ps
   massSamp[ii]=rob_p_mean
   massErrSamp[ii]=rob_p_sigma

   printf,u,ii,nDiffSamp[ii],zSamp[ii],LXSamp[ii],mLXSamp[ii],massSamp[ii],massErrSamp[ii]
endfor

close,u
free_lun,u

stop
end
