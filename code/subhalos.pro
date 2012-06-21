function zMatchField, zMem, zField, zMin, zMax
; select a z-matched field sample to the member sample
; returns array of indices sampled from zField array

binWidth=0.02

yMem=histogram(zMem,min=zMin,max=zMax,binsize=binWidth)
yField=histogram(zField,min=zMin,max=zMax,binsize=binWidth,reverse_indices=ri)

;; to ensure we have enough field galaxies in each bin, normalize
;; distribution where Nfield/Nmem is lowest
;norm=min(float(yField)/yMem) ; should return finite value

; normalize to get an equal number of field galaxies as members
norm=1.

if(norm LT 1) then begin
   print,'problem with z-matched field sample; there should be more field galaxies than members at all redshifts'
   stop
endif

nSubField=yMem*norm ; this is the number of field galaxies in each z-bin we should end up with
for ii=0,n_elements(yMem)-1 do begin
   ; first get the indices of field members in this zBin, using histogram magic
   ; see http://www.astro.virginia.edu/class/oconnell/astr511/IDLresources/histogram_tutorial-JDSmith.html
   if(ri[ii] NE ri[ii+1]) then tmpInd=ri[ri[ii]:ri[ii+1]-1]
   ; now pick a random subset of these indices, with the number selected                     ; scaled to nSubField
   ; see http://www.idlcoyote.com/code_tips/randomindex.html
   if(nSubField[ii] GT 0) then sel=(sort(randomu(seed,n_elements(tmpInd))))[0:nSubField[ii]-1]

   if(n_elements(selField) EQ 0 AND n_elements(sel) GT 0) then selField=tmpInd[sel] $
   else if(n_elements(sel) GT 0) then selField=[selField,tmpInd[sel]] ; kludge to accumulate array
   undefine,sel ; clear sel so it doesn't get used next time
endfor

return, selField
end

pro add_group_centers,acs
  alpha_mmgg_scale_specz=fltarr(n_elements(acs))
  delta_mmgg_scale_specz=fltarr(n_elements(acs))
  for ii=0,n_elements(acs)-1 do begin
     cenInd=where(acs.mmgg_scale_specz EQ 1 $
                  AND acs.group_id_best_specz EQ acs[ii].group_id_best_specz $
                  ,nCen)
     if(nCen GT 0) then begin
        alpha_mmgg_scale_specz[ii]=acs[cenInd].alpha_j2000
        delta_mmgg_scale_specz[ii]=acs[cenInd].delta_j2000
     endif else begin
        alpha_mmgg_scale_specz[ii]=-999.
        delta_mmgg_scale_specz[ii]=-999.
     endelse
  endfor
  acs=mrg_addcol(acs,'ALPHA_MMGG_SCALE_SPECZ',alpha_mmgg_scale_specz)
  acs=mrg_addcol(acs,'DELTA_MMGG_SCALE_SPECZ',delta_mmgg_scale_specz)
end

pro create_subhalo_cats, catFiles, acs

inner=where(acs.p_mem_best_specz GT 0.5 $
            AND acs.group_flag_best_specz EQ 1 $
            AND acs.mmgg_scale_specz EQ 0 $
            AND acs.mmgg_scale EQ 0 $
            AND acs.dist_bcg_r200_specz LE 0.45 $
            ,nInner)

outer=where(acs.p_mem_best_specz GT 0.5 $
            AND acs.group_flag_best_specz EQ 1 $
            AND acs.mmgg_scale_specz EQ 0 $
            AND acs.mmgg_scale EQ 0 $
            AND acs.dist_bcg_r200_specz GT 0.45 $
            ,nOuter)

field=where(acs.p_mem_best_specz LE 0 $
            AND acs.mmgg_scale_specz EQ 0 $
            AND acs.mmgg_scale EQ 0 $
            AND acs.mag_auto LT 24.2 $
            AND acs.zphot GT 0 $
            AND acs.zphot LT 1 $
            AND acs.group_projection_mmgg_specz EQ 0 $
            ,nField)

if(nOuter GT nInner) then begin
;   subsel=zMatchField(acs[inner].zphot,acs[outer].zphot,0.,1.)
   subsel=(sort(randomu(seed,nOuter)))[0:nInner-1]
   outer=outer[subsel]
endif else if(nInner GT nOuter) then begin
;   subsel=zMatchField(acs[outer].zphot,acs[inner].zphot,0.,1.)
   subsel=(sort(randomu(seed,nInner)))[0:nOuter-1]
   inner=inner[subsel]
endif
subsel=zMatchField(acs[inner].zphot,acs[field].zphot,0.,1.)
field=field[subsel] ; field shoud now have same number of elements as inner

print,'inner, outer, field'
print,'N',n_elements(inner),n_elements(outer),n_elements(field)
print,'mean M*',mean(acs[inner].kevin_mstar),mean(acs[outer].kevin_mstar),mean(acs[field].kevin_mstar)
print,'median M*',median(acs[inner].kevin_mstar),median(acs[outer].kevin_mstar),median(acs[field].kevin_mstar)
print,'stddev M*',stddev(acs[inner].kevin_mstar),stddev(acs[outer].kevin_mstar),stddev(acs[field].kevin_mstar)
print,'mean z',mean(acs[inner].zphot),mean(acs[outer].zphot),mean(acs[field].zphot)

acs=mrg_addcol(acs,'LX_SCALE',replicate(0.,n_elements(acs)))
acs=mrg_addcol(acs,'LENSING_M200',replicate(0.,n_elements(acs)))
acs=mrg_addcol(acs,'LENSING_R200_MPC',replicate(0.,n_elements(acs)))
acs=mrg_addcol(acs,'BOX',replicate(0.,n_elements(acs)))
acs=mrg_addcol(acs,'ID',replicate(0.,n_elements(acs)))
acs=mrg_addcol(acs,'FLAG_INCLUDE',replicate(1.,n_elements(acs)))

add_group_centers, acs

mwrfits,acs[inner],catFiles[0],/create
mwrfits,acs[outer],catFiles[1],/create
mwrfits,acs[field],catFiles[2],/create

for ii=0, n_elements(catFiles)-1 do begin
   assign_boxes,catFiles[ii],catFiles[ii]
endfor

return
end

pro subhalos

acs=mrdfits("~/data/cosmos/code/lensing18_20110914.fits",1,$
           columns=['ALPHA_J2000',$
                    'DELTA_J2000',$
                    'IDENT',$
                    'ZPHOT',$
                    'MAG_AUTO',$
                    'KEVIN_MSTAR',$
                    'P_MEM_BEST_SPECZ',$
                    'GROUP_ID_BEST_SPECZ',$
                    'MMGG_SCALE_SPECZ',$
                    'MMGG_SCALE',$
                    'GROUP_FLAG_BEST_SPECZ',$
                    'DIST_BCG_R200_SPECZ',$
                    'GROUP_PROJECTION_MMGG_SPECZ'$
                   ])
group=mrdfits("~/data/cosmos/code/group5_20110914.fits",1)

outDir="~/data/cosmos/groups_lensing/outfiles/subhalos/"
plotDir="~/data/cosmos/groups_lensing/plots/subhalos/"
catNames=["inner","outer","field"]
catFiles=outDir+catNames+"Cat.fits"
lensOutFiles=outDir+"lens_"+catNames+".fits"
lensPlotFiles=plotDir+catNames+".eps"
lensPlotModelFiles=plotDir+catNames+"_models.eps"

sel=where(acs.zphot GT 0 $
          AND acs.zphot LE 1 $
          AND acs.mag_auto LT 24.2 $
          AND acs.kevin_mstar GT 10.3 $
          AND acs.kevin_mstar LT 11.5 $
         )
acs=acs[sel]

; flag group members where centers are unclear
for ii=0,n_elements(group)-1 do begin
   mem=where(acs.group_id_best_specz EQ group[ii].id, nMem)
   if(nMem GT 0) then begin
      if(group[ii].id_mmgg_scale_specz NE group[ii].id_mmgg_scale $
         OR group[ii].id_mmgg_scale_specz NE group[ii].id_mmgg_r200 $
         OR group[ii].id_mmgg_scale_specz NE group[ii].id_mlgg_r200 $
         OR group[ii].id_mmgg_scale_specz NE group[ii].id_mlgg_scale) then $
            acs[mem].group_flag_best_specz=0
   endif
endfor

create_subhalo_cats,catFiles, acs

infile_source='/Users/alexie/Work/Weak_lensing/GG_cat_2006/gglensing_source_v1.7.fits' ; Using the new catalog (photoz version 1.7)

innerRadiusKpc=20.
secondRadiusKpc=70.
maxRadiusKpc=1000.
nRadiusBins=7.
minLensZ=0.
maxLensZ=1.
minLensMass=0.
maxLensMass=15.
box_factor=20
zscheme=2
center='galaxy'
refcen='mmgg_scale_specz'

cen_type='ps'

for ii=0,n_elements(catNames)-1 do begin
   run_gg_offset, infile_source, catFiles[ii], lensOutFiles[ii], innerRadiusKpc, secondRadiusKpc, maxRadiusKpc, nRadiusBins, minLensZ, maxLensZ, minLensMass, maxLensMass, box_factor, zscheme, /xgroups,/usespecz,center=center,/emp_var

   if(ii EQ 0 OR ii EQ 1) then begin
      fit_t = [$
              2,$               ; 0  M0    : baryonic mass
              1,$               ; 1  R_vir : NFW virial mass
              2,$               ; 2  C     : NFW concentration
              0,$               ; 3  alpha : fraction
              0,$               ; 4  bias
              0,$               ; 5  m_sigma
              0 ]               ; 6  offset

      p_mean=[13.5]

      plot_lensing_results, lensOutFiles[ii], lensPlotFiles[ii], p_mean, fit_t, cen_type, $
                            center=center,refcen=refcen,groupFile=catFiles[ii],/use_m200
      plot_lensing_results, lensOutFiles[ii], lensPlotModelFiles[ii], p_mean, fit_t, cen_type, $
                            center=center,refcen=refcen,groupFile=catFiles[ii],/use_m200,/models
   endif else begin
      fit_t = [$
              2,$               ; 0  M0    : baryonic mass
              0,$               ; 1  R_vir : NFW virial mass
              0,$               ; 2  C     : NFW concentration
              0,$               ; 3  alpha : fraction
              0,$               ; 4  bias
              0,$               ; 5  m_sigma
              0 ]               ; 6  offset
      p_mean=-1
      plot_lensing_results, lensOutFiles[ii], lensPlotFiles[ii], p_mean, fit_t, cen_type, $
                            /use_m200
      plot_lensing_results, lensOutFiles[ii], lensPlotModelFiles[ii], p_mean, fit_t, cen_type, $
                            /use_m200,/models
   endelse
endfor

end
