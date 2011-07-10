pro test_al_mass,nostackx=nostackx,noemp=noemp

; try to replicate masses from Table 1 of AL's LX-M paper. 
; motivated by results from jackknife sampling which tended to show
; lower mean WL masses than the LX relation predicts.

; first define similar LX-z bins as in Fig. 3 and write out catalog
; files then do the stacked lensing

; use old catalog to keep ALs centers
;gOld=mrdfits("~/data/cosmos/catalogs/group11_0126.fits",1)
gOld=mrdfits("~alexie/Work/Weak_lensing/Group_cat_2008/group11.fits",1)

infile_source='/Users/alexie/Work/Weak_lensing/GG_cat_2006/gglensing_source_v1.7.fits' ; Using the new catalog (photoz version 1.7)

fileExt=''
if(keyword_set(nostackx) AND keyword_set(noemp)) then $
   fileExt='_nostackx_noemp' $
else if (keyword_set(nostackx)) then $
   fileExt='_nostackx' $
else if (keyword_set(noemp)) then $
   fileExt='_noemp'


; Set paths for output files
outDir='/Users/mgeorge/data/cosmos/groups_lensing/outfiles/test_al_mass/'
if(NOT(file_test(outDir))) then file_mkdir,outDir
outTableFile=outDir+'table1'+fileExt+'.txt'
openw,u,outTableFile,/get_lun,width=1000
printf,u,'# BinID  Nlens  LX_scale  Mlx  Mwl  Mwl_err  z'
printf,u,'#  #       #    x10^42       x10^13         '

;AL's bins
; cf. get_group_paper_boxes()
lxBins=[ $
       [43.30,43.70],$          ;A0
       [43.01,43.41],$          ;A1
       [42.61,43.01],$          ;A2
       [42.21,42.61],$          ;A3
       [41.71,42.11],$          ;A4
       [41.81,42.21],$          ;A5
       [42.30,42.90],$          ;A6
       [42.37,42.97],$          ;A7
       [42.75,43.35],$          ;A8
       [43.30,43.70],$ ;A9 ; repeat the first bin, we can use this to test the repeatability of the lensing code when a different number of identical groups are used
       [43.30,43.70],$ ;A10
       [43.30,43.70],$ ;A11
       [43.30,43.70],$ ;A12
       [43.30,43.70],$ ;A13
       [43.30,43.70],$ ;A14
       [43.30,43.70],$ ;A15
       [43.30,43.70] $ ;A16
       ]
zBins=[ $
      [0.17,0.27],$                            ;A0
      [0.3,0.4],$                              ;A1
      [0.3,0.4],$                              ;A2
      [0.3,0.4],$                              ;A3
      [0.18,0.28],$                            ;A4
      [0.3,0.4],$                              ;A5
      [0.43,0.58],$                            ;A6
      [0.595,0.745],$                          ;A7
      [0.816,0.966],$                          ;A8
      [0.17,0.27],$                            ;A9
      [0.17,0.27],$                            ;A10
      [0.17,0.27],$                            ;A11
      [0.17,0.27],$                            ;A12
      [0.17,0.27],$                            ;A13
      [0.17,0.27],$                            ;A14
      [0.17,0.27],$                            ;A15
      [0.17,0.27] $                            ;A16
      ]
nRep=[2,1,1,1,1,1,1,1,1,2,2,2,10,20,30,50,100] ; number of times to repeat single group for bins A0,A9+
nBins=(size(lxBins,/dim))[1]

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
minRadiusKpc=50.
maxRadiusKpc=2000.
nRadiusBins=8

cenName='bcg1'

; arrays to record mean and stddev mass from mcmc
massMean=fltarr(nBins)
massErr=fltarr(nBins)

for ii=0,nBins-1 do begin

   ; select groups and write catalog files
   sel=where(gOld.zphot GE zBins[0,ii] $
             AND gOld.zphot LE zBins[1,ii] $
             AND gOld.lx_scale GE lxBins[0,ii] $
             AND gOld.lx_scale LE lxBins[1,ii] $
             AND gOld.bcg_flag GE 2 $ ; AL used BCG_FLAG=2 or 3 in her table 1
             ,nSel)
   if(nSel EQ 1) then begin
      sel=replicate(sel,nRep[ii])      ; duplicate entry to avoid scalar/array bugs
      nSel=n_elements(sel)
   endif
   if(nSel EQ 0) then begin
      print,'no groups in bin',ii
      stop
   endif

   ; add tags missing from old group file
   ; that are needed in new lensing code
   gNew=create_struct(gOld[0],$
                      'FLAG_INCLUDE',1)
   gNew=replicate(gNew,nSel)

   for jj=0,n_tags(gOld)-1 do begin ; fill in each tag jj with the old info
                                ; note - we don't need to know the tag names
      gNew.(jj)=gOld[sel].(jj)
   endfor
   
   ; write out new group file for this bin
   groupFile=outDir+'group_a'+string(ii,format='(I2)')+fileExt+'.fits'
   mwrfits,gNew,groupFile,/create


   ; LENSING
   lensOutFile=strcompress(outDir+'center_a'+string(ii,format='(I2)')+fileExt+'.fits',/remove_all)
   run_gg_offset, infile_source, groupFile, lensOutFile, minRadiusKpc, maxRadiusKpc, nRadiusBins, minLensZ, maxLensZ, minLensMass, maxLensMass, box_factor, zscheme, /xgroups,/usespecz,center=cenName,stackx=~(keyword_set(nostackx)),emp_var=~(keyword_set(noemp))

   ; Fit model to the lensing signal
   run_ds_mcmc, lensOutFile, fit_t, rob_p_mean, rob_p_sigma, stackx=~(keyword_set(nostackx))
   massMean[ii]=rob_p_mean
   massErr[ii]=rob_p_sigma

   
   printf,u,'a'+string(ii,format='(I2)'),nSel,mean(10.^(gnew.lx_scale-42.)),mean(10.^(gnew.lensing_m200-13.)),10.^(massMean[ii]-13.),10.^(massErr[ii]),mean(gnew.zphot)

endfor
close,u
free_lun,u

stop
end
