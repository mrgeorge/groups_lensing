function nfw_ds_offset, r, p, zl, groupFile, $
                        r200 = r200, $
                        r180 = r180, $
                        center=center, $
                        refcen=refcen, $
                        roff=roff

; center is the string name of the center under consideration
; refcen is the best halo center, default 'mmgg_scale'
; roff is the mean offset in Mpc - this option can be set to specify
;                                  the offset explicitly. If not set,
;                                  the offset distribution is taken to be the
;                                  distribution between center and refcen
; r is distance in MPC
; p[0] is r200 the virial (or 200) radius in MPC
; p[1] is c the concentration parameter
; r180   -> 180 x mean density
; r200   -> 200 x critical density

if(n_elements(roff) GT 0) then begin
   if(min(roff) LT 0.) then begin
      print,'nfw_ds_offset returning ds=0: roff=',roff
      return, 0.*r              ; return 0's since roff shouldn't be negative
   endif
   if(max(roff) EQ 0.) then return, nfw_ds(r,p,zl,r200=r200,r180=r180) ; no offsets

   ; (else)
   ngroups=n_elements(roff)
endif else begin
   if(NOT(keyword_set(refcen))) then refcen='mmgg_scale'

   ; read in group file to get offsets
   group=mrdfits(groupFile,1)

   ; get array of offset lengths
   groupsel=group[where(group.flag_include EQ 1)]
   get_center_coords,groupsel,center,ra1,dec1,sel1
   get_center_coords,groupsel,refcen,ra2,dec2,sel2
   roff=distance(ra1,dec1,ra2,dec2)*3600.*groupsel.lensing_r200_mpc/groupsel.lensing_r200_as ; Mpc
   offset_kpc=roff*1000.
   min_offset_kpc=50.
   sel=where(ra1 GT 0 $
             AND ra2 GT 0 $
             AND offset_kpc GT min_offset_kpc $
            )
   roff=roff[sel]
   ngroups=n_elements(roff)
endelse

; get relative weights from number of sources in bin behind lens
; (need to modify gglensing_v6 to get this - use 1 for now)
groupWeight=replicate(1./ngroups,ngroups)

; loop over groups
sigmaR=dblarr(n_elements(r),ngroups)
sigmaMean=dblarr(n_elements(r),ngroups)
deltaSigma=dblarr(n_elements(r))
for g=0,ngroups-1 do begin
   ; call AL's version
   deltaSigma+=groupWeight[g]*nfw_ds_offset_sat(r,p,zl,roff[g],r200=r200,r180=r180)
endfor

return, deltaSigma
end
