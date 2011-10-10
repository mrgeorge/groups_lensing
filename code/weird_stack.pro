pro do_weird_stack, groupSel, cenName, selType
; write a group file
; measure the lensing signal
; plot the lensing signal

; groupSel is group[sel], i.e. a subarray of group structures
; cenName is the name of the center used, e.g. 'mlgg_r200'
; selType is a description of the selection of groups, for filenames

; Set paths for output files
dirName='weird_stack'
fileDir='~/data/cosmos/groups_lensing/outfiles/'+dirName+'/bin_10_70_1000_7/'
plotDir='~/data/cosmos/groups_lensing/plots/'+dirName+'/bin_10_70_1000_7/'
if(NOT(file_test(fileDir))) then file_mkdir,fileDir
if(NOT(file_test(plotDir))) then file_mkdir,plotDir

groupFile=fileDir+'group_'+selType+'.fits'
lensOutFile=fileDir+'center_'+selType+'.fits'
infile_source='/Users/alexie/Work/Weak_lensing/GG_cat_2006/gglensing_source_v1.7.fits' ; Using the new catalog (photoz version 1.7)

mwrfits,groupSel,groupFile,/create

innerRadiusKpc=10.
secondRadiusKpc=70.
maxRadiusKpc=1000.
nRadiusBins=7
minLensZ=0.
maxLensZ=1.0
minLensMass=12.
maxLensMass=15.
box_factor=20.
zscheme=2

run_gg_offset, infile_source, groupFile, lensOutFile, innerRadiusKpc, secondRadiusKpc, maxRadiusKpc, nRadiusBins, minLensZ, maxLensZ, minLensMass, maxLensMass, box_factor,zscheme,/xgroups,/usespecz,center=cenName ; no stackx or emp_var

fitType = [$
2,$             ; 0  M0    : baryonic mass
1,$             ; 1  R_vir : NFW virial mass
2,$             ; 2  C     : NFW concentration
0,$             ; 5  alpha : fraction
0,$             ; 6  bias
0 ]             ; 6  m_sigma
run_ds_mcmc, lensOutFile, fitType, rob_p_mean, rob_p_sigma
plot_lensing_results,lensOutFile,plotDir+selType,rob_p_mean,fitType,/use_m200,/models
end

pro weird_stack

; plot lensing signal for stacks of strange groups
; e.g., blue centers, cases where mass gap is smaller than uncertainty,
;       different combinations of of centers where they disagree

groupAll=mrdfits("~alexie/Work/GroupCatalogs/cosmos_xgroups_20110209.fits",1)
sel=where(groupAll.flag_include EQ 1,nGroups)
group=groupAll[sel]

acs=mrdfits("~alexie/Work/GroupCatalogs/lensing15_20110209.fits",1)

mmggScaleInd=lonarr(nGroups)
mmgg200Ind=lonarr(nGroups)
mlggScaleInd=lonarr(nGroups)
mlgg200Ind=lonarr(nGroups)

for ii=0,nGroups-1 do begin
   mmggScaleInd[ii]=where(acs.ident EQ group[ii].id_mmgg_scale)
   mmgg200Ind[ii]=where(acs.ident EQ group[ii].id_mmgg_r200)
   mlggScaleInd[ii]=where(acs.ident EQ group[ii].id_mlgg_scale)
   mlgg200Ind[ii]=where(acs.ident EQ group[ii].id_mlgg_r200)
endfor

; define centers by color
blueMLGG200=where(acs[mlgg200Ind].mnuv_mr LT 1.2) ; 21
greenMLGG200=where(acs[mlgg200Ind].mnuv_mr GT 1.2 AND acs[mlgg200Ind].mnuv_mr LT 3.5) ; 22
redMLGG200=where(acs[mlgg200Ind].mnuv_mr GT 3.5) ; 86
blueMMGG200=where(acs[mmgg200Ind].mnuv_mr LT 1.2) ; 8
greenMMGG200=where(acs[mmgg200Ind].mnuv_mr GT 1.2 AND acs[mmgg200Ind].mnuv_mr LT 3.5) ; 24
redMMGG200=where(acs[mmgg200Ind].mnuv_mr GT 3.5) ; 97
blueMMGGscale=where(acs[mmggScaleInd].mnuv_mr LT 1.2) ; 8
greenMMGGscale=where(acs[mmggScaleInd].mnuv_mr GT 1.2 AND acs[mmggScaleInd].mnuv_mr LT 3.5) ; 24
redMMGGscale=where(acs[mmggScaleInd].mnuv_mr GT 3.5) ; 97

;do_weird_stack,group[blueMLGG200],'mlgg_r200','blue_mlgg200'
;do_weird_stack,group[greenMLGG200],'mlgg_r200','green_mlgg200'
;do_weird_stack,group[redMLGG200],'mlgg_r200','red_mlgg200'
;do_weird_stack,group[blueMMGG200],'mmgg_r200','blue_mmgg200'
;do_weird_stack,group[greenMMGG200],'mmgg_r200','green_mmgg200'
;do_weird_stack,group[redMMGG200],'mmgg_r200','red_mmgg200'
;do_weird_stack,group[blueMMGGscale],'mmgg_scale','blue_mmgg_scale'
;do_weird_stack,group[greenMMGGscale],'mmgg_scale','green_mmgg_scale'
;do_weird_stack,group[redMMGGscale],'mmgg_scale','red_mmgg_scale'

; define centers by morphology
; note "late" includes lates, irregulars, and a few failed or "other" classifications
earlyMLGG200=where(acs[mlgg200Ind].zest_type EQ 1 OR (acs[mlgg200Ind].zest_type EQ 2 AND acs[mlgg200Ind].zest_bulge EQ 0),complement=lateMLGG200)
earlyMMGG200=where(acs[mmgg200Ind].zest_type EQ 1 OR (acs[mmgg200Ind].zest_type EQ 2 AND acs[mmgg200Ind].zest_bulge EQ 0),complement=lateMMGG200)
earlyMMGGscale=where(acs[mmggScaleInd].zest_type EQ 1 OR (acs[mmggScaleInd].zest_type EQ 2 AND acs[mmggScaleInd].zest_bulge EQ 0),complement=lateMMGGscale)

;do_weird_stack,group[earlyMLGG200],'mlgg_r200','early_mlgg200'
;do_weird_stack,group[lateMLGG200],'mlgg_r200','late_mlgg200'
;do_weird_stack,group[earlyMMGG200],'mmgg_r200','early_mmgg200'
;do_weird_stack,group[lateMMGG200],'mmgg_r200','late_mmgg200'
;do_weird_stack,group[earlyMMGGscale],'mmgg_scale','early_mmgg_scale'
;do_weird_stack,group[lateMMGGscale],'mmgg_scale','late_mmgg_scale'

; select groups with a nearby group in projection
dis_as=fltarr(nGroups)
for ii=0,nGroups-1 do begin
   others=where(groupAll.id NE group[ii].id)
   dis_as[ii]=min(distance(groupAll[others].alpha_ellipse,groupAll[others].delta_ellipse,group[ii].alpha_mmgg_scale,group[ii].delta_mmgg_scale)*3600)
endfor
nearby=where(dis_as LT 1.5*group.lensing_r200_as,complement=far)

;do_weird_stack,group[nearby],'mmgg_scale','nearby_mmgg_scale'
;do_weird_stack,group[far],'mmgg_scale','far_mmgg_scale'


;;; FOCUSING ON GROUPS WITH DISCREPANT GALAXY CENTERS

; define visual center from inspection notes
readcol,'mrg_center_notes_20110808_distilled.txt',id,imnum,flag,msm2,msls,msl2,m2ls,m2l2,lsl2,best,conf,merg,proj,bad,format='i,i,i,i,i,i,i,i,i,a,i,i,i,i'
match,id,group.id,m1,m2
groupWeird=group[m2] 
groupNonWeird=group[where(group.id_mmgg_scale EQ group.id_mmgg_r200 $
                          AND group.id_mmgg_scale EQ group.id_mlgg_scale $
                          AND group.id_mmgg_scale EQ group.id_mlgg_r200 $
                          AND group.id_mmgg_r200 EQ group.id_mlgg_scale $
                          AND group.id_mmgg_r200 EQ group.id_mlgg_r200 $
                          AND group.id_mlgg_scale EQ group.id_mlgg_r200)]

alpha_visual=fltarr(n_elements(m1))
delta_visual=fltarr(n_elements(m1))
id_visual=fltarr(n_elements(m1))
visual_mstar=fltarr(n_elements(m1))
for ii=0,n_elements(m1)-1 do begin
   case best[m1[ii]] of
      'Ms': begin
         alpha_visual[ii]=groupWeird[ii].alpha_mmgg_scale
         delta_visual[ii]=groupWeird[ii].delta_mmgg_scale
         id_visual[ii]=groupWeird[ii].id_mmgg_scale
         visual_mstar[ii]=groupWeird[ii].mmgg_scale_mstar
      end
      'M2': begin
         alpha_visual[ii]=groupWeird[ii].alpha_mmgg_r200
         delta_visual[ii]=groupWeird[ii].delta_mmgg_r200
         id_visual[ii]=groupWeird[ii].id_mmgg_r200
         visual_mstar[ii]=groupWeird[ii].mmgg_r200_mstar
      end
      'Ls': begin
         alpha_visual[ii]=groupWeird[ii].alpha_mlgg_scale
         delta_visual[ii]=groupWeird[ii].delta_mlgg_scale
         id_visual[ii]=groupWeird[ii].id_mlgg_scale
         visual_mstar[ii]=groupWeird[ii].mlgg_scale_mstar
      end
      'L2': begin
         alpha_visual[ii]=groupWeird[ii].alpha_mlgg_r200
         delta_visual[ii]=groupWeird[ii].delta_mlgg_r200
         id_visual[ii]=groupWeird[ii].id_mlgg_r200
         visual_mstar[ii]=groupWeird[ii].mlgg_r200_mstar
      end
      else: begin
         galID=long(best[m1[ii]])
         if(galID LE 0) then stop ; best identifier is unknown
         galInd=where(acs.ident EQ galID)
         alpha_visual[ii]=acs[galInd].alpha_j2000
         delta_visual[ii]=acs[galInd].delta_j2000
         id_visual[ii]=galID
         visual_mstar[ii]=acs[galInd].kevin_mstar
      endelse
   endcase
endfor
groupWeird=mrg_addcol(groupWeird,'ALPHA_VISUAL',alpha_visual)
groupWeird=mrg_addcol(groupWeird,'DELTA_VISUAL',delta_visual)
groupWeird=mrg_addcol(groupWeird,'ID_VISUAL',id_visual)
groupWeird=mrg_addcol(groupWeird,'VISUAL_MSTAR',visual_mstar)

; Stack on different centers for this sample of groups which has
; discrepant centers that were visually inspected
do_weird_stack,groupWeird,'visual','weird_visual'
do_weird_stack,groupWeird,'mmgg_scale','weird_mmgg_scale'
do_weird_stack,groupWeird,'mlgg_scale','weird_mlgg_scale'
do_weird_stack,groupWeird,'mmgg_r200','weird_mmgg_r200'
do_weird_stack,groupWeird,'mlgg_r200','weird_mlgg_r200'
do_weird_stack,groupWeird,'cn','weird_cn'
do_weird_stack,groupWeird,'cm','weird_cm'
do_weird_stack,groupWeird,'cl','weird_cl'

; try getting rid of different cases like mergers or possible projections
nonBad=where(bad[m1] EQ 0)
nonProj=where(proj[m1] EQ 0)
nonMerg=where(merg[m1] EQ 0)
confident=where(conf[m1] LE 2)
do_weird_stack,groupWeird[nonBad],'visual','weird_nonBad_visual'
do_weird_stack,groupWeird[nonProj],'visual','weird_nonProj_visual'
do_weird_stack,groupWeird[nonMerg],'visual','weird_nonMerg_visual'
do_weird_stack,groupWeird[confident],'visual','weird_confident_visual'


; compare with stack of groups where all gal centers agree
do_weird_stack,groupNonWeird,'mmgg_scale','nonWeird_mmgg_scale'
;do_weird_stack,groupNonWeird,'mmgg_r200','nonWeird_mmgg_r200'
;do_weird_stack,groupNonWeird,'mlgg_scale','nonWeird_mlgg_scale'
;do_weird_stack,groupNonWeird,'mlgg_r200','nonWeird_mlgg_r200'

; try to get a cleaner signal from the nonWeird sample
hiM=where(groupNonWeird.mmgg_scale_mstar GT 11.)
do_weird_stack,groupNonWeird[hiM],'mmgg_scale','nonWeirdHiM_mmgg_scale'

gap=fltarr(n_elements(groupNonWeird))
for ii=0,n_elements(gap)-1 do begin
   mem=where(acs.p_mem_best GT 0.5 AND acs.group_id_best EQ groupNonWeird[ii].id)
   topTwo=acs[mem[(reverse(sort(acs[mem].kevin_mstar)))[0:1]]].kevin_mstar
   gap[ii]=topTwo[0]-topTwo[1]
endfor   
bigGap=where(gap GT 0.3,complement=smallGap)

do_weird_stack,groupNonWeird[bigGap],'mmgg_scale','nonWeirdBigGap_mmgg_scale'
do_weird_stack,groupNonWeird[smallGap],'mmgg_scale','nonWeirdSmallGap_mmgg_scale'

rich=where(groupNonWeird.n_mem LT 20,complement=poor)

do_weird_stack,groupNonWeird[rich],'mmgg_scale','nonWeirdRich_mmgg_scale'
do_weird_stack,groupNonWeird[poor],'mmgg_scale','nonWeirdPoor_mmgg_scale'


;loz=where(groupNonWeird.zphot GT 0.2 $
;          AND groupNonWeird.zphot LT 0.6)
;hiz=where(groupNonWeird.zphot GT 0.6 $
;          AND groupNonWeird.zphot LT 1.0)
loz=where(group.zphot GT 0.2 $
          AND group.zphot LT 0.6)
hiz=where(group.zphot GT 0.6 $
          AND group.zphot LT 1.0)
allz=where(group.zphot GT 0.2 $
           AND group.zphot LT 1.0)
;do_weird_stack,group[loz],'mmgg_scale','loZ_mmgg_scale'
;do_weird_stack,group[hiz],'mmgg_scale','hiZ_mmgg_scale'
;do_weird_stack,group[allz],'mmgg_scale','allZ_mmgg_scale'


; compare centers for groups where X-ray and MMGG_scale are offset
; more than the X-ray uncertainty (and more than 50 kpc)
d_x_ms_kpc=distance(group.alpha_ellipse,group.delta_ellipse,group.alpha_mmgg_scale,group.delta_mmgg_scale)*3600.*group.lensing_r200_mpc*1000./group.lensing_r200_as
xoff=where(group.flag1 EQ 1 $
           AND d_x_ms_kpc GT 50. $
           AND group.pos_err_ellipse*3600.*group.lensing_r200_mpc*1000./group.lensing_r200_as LT 50.)
do_weird_stack,group[xoff],'mmgg_scale','xoff_mmgg_scale'
do_weird_stack,group[xoff],'xray','xoff_xray'

; repeat using AF's centers (alpha_j2000) instead of ellipse
d_xaf_ms_kpc=distance(group.alpha_j2000,group.delta_j2000,group.alpha_mmgg_scale,group.delta_mmgg_scale)*3600.*group.lensing_r200_mpc*1000./group.lensing_r200_as
xaf_off=where(group.flag1 EQ 1 $
           AND d_xaf_ms_kpc GT 50. $
           AND group.pos_err_ellipse*3600.*group.lensing_r200_mpc*1000./group.lensing_r200_as LT 50.)
do_weird_stack,group[xaf_off],'mmgg_scale','xaf_off_mmgg_scale'
do_weird_stack,group[xaf_off],'af','xaf_off_af'

; AL saw problems with a hi-M sample
hiM=where(group.zphot GT 0.2 $
          AND group.zphot LT 1. $
          AND group.lensing_m200 GT 13.55)
do_weird_stack,group[hiM],'mmgg_scale','hiM_mmgg_scale'

hiMClean=where(group.zphot GT 0.2 $
               AND group.zphot LT 1. $
               AND group.lensing_m200 GT 13.55 $
               AND group.id_mmgg_scale EQ group.id_mlgg_scale)
do_weird_stack,group[hiMClean],'mmgg_scale','hiMClean_mmgg_scale'

end
