pro jackknife_stack

; stack random subsamples of groups to determine the distribution of
; masses obtained. This is a test to see whether groups with
; MMGG_scale != MMGG_r200 have anomalously low masses.

group=mrdfits("~alexie/Work/GroupCatalogs/cosmos_xgroups_20110209.fits",1)
sel=where(group.flag_include EQ 1)
group=group[sel]

bad=where(group.id_mmgg_scale NE group.id_mmgg_r200)
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

minRadiusKpc=50.
maxRadiusKpc=2000.
nRadiusBins=8
stackx=0
emp_var=0
cenName='mmgg_scale'
dirName='bin_'+string(minRadiusKpc,format='(I0)')+'_'+string(maxRadiusKpc,format='(I0)')+'_'+string(nRadiusBins,format='(I0)')
if(keyword_set(stackx)) then dirName += '_sx'
if(keyword_set(emp_var)) then dirName += '_emp'
fileDir='~/data/cosmos/groups_lensing/outfiles/'+dirName+'/'
lensOutFile=strcompress(fileDir+'center_'+cenName+'+_jackknife.fits',/remove_all)
sampGroupFile=fileDir+'group_jackknife.fits'
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
0,$             ; 5  alpha : fraction
0,$             ; 6  bias
0 ]             ; 6  m_sigma

openw,u,fileDir+'jackknife_results.txt',/get_lun,width=1000
printf,u,'# Jackknife sample size',sampleSize
printf,u,'# Run  nDiff  z  LX  mLX  mWL  sig_mWL'

for ii=0,nSamples-1 do begin
   print,'==========================='
   print,'Sample ',ii,' of ',nSamples
   print,'==========================='

   ; choose N random indices from Ngroups without replacement 
   rand=randomu(seed,n_elements(group))
   sel=(sort(rand))[0:sampleSize-1]
   mwrfits,group[sel],sampGroupFile,/create 

   LXSamp[ii]=median(group[sel].lx_scale)
   mLXSamp[ii]=median(group[sel].lensing_m200)
   zSamp[ii]=median(group[sel].zphot)
   nDiffSamp[ii]=n_elements(where(group[sel].id_mmgg_scale NE group[sel].id_mmgg_r200))

   run_gg_offset, infile_source, sampGroupFile, lensOutFile, minRadiusKpc, maxRadiusKpc, nRadiusBins, minLensZ, maxLensZ, minLensMass, maxLensMass, box_factor, zscheme, /xgroups,/usespecz,center=cenName,stackx=keyword_set(stackx),emp_var=keyword_set(emp_var)

   ; Fit model to the lensing signal
   run_ds_mcmc, lensOutFile, fitType, rob_p_mean, rob_p_sigma
   massSamp[ii]=rob_p_mean
   massErrSamp[ii]=rob_p_sigma

   printf,u,ii,nDiffSamp[ii],zSamp[ii],LXSamp[ii],mLXSamp[ii],massSamp[ii],massErrSamp[ii]
endfor

close,u
free_lun,u

stop
end
