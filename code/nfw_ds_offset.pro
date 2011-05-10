function nfw_ds_offset, r, p, zl, $
                        r200 = r200, $
                        r180 = r180, $
                        center=center, $
                        refcen=refcen 

; center is the string name of the center under consideration
; refcen is the best halo center, default 'mmgg_scale'
; r is distance in MPC
; p[0] is r200 the virial (or 200) radius in MPC
; p[1] is c the concentration parameter
; r180   -> 180 x mean density
; r200   -> 200 x critical density

if(NOT(keyword_set(refcen))) then refcen='mmgg_scale'

common xgroup_cat, group
if(NOT(keyword_set(group))) then begin
   groupfile="~alexie/Work/Weak_lensing/Group_cat_june_2010/group6.fits"
   group=mrdfits(groupfile,1)
endif

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

; get relative weights from number of sources in bin behind lens
; (need to modify gglensing_v6 to get this - use 1 for now)
groupWeight=replicate(1./ngroups,ngroups)

; loop over groups
sigmaR=dblarr(n_elements(r),ngroups)
sigmaMean=dblarr(n_elements(r),ngroups)
for g=0,ngroups-1 do begin

   ; \Sigma(R|R_{off})
   sigmaR[*,g] = groupWeight[g] * nfw_sigma_offset(r,roff[g],p,zl,r200=keyword_set(r200),r180=keyword_set(r180))

   ; \bar{\Sigma}(<R|R_{off})
   for i=0,n_elements(r)-1 do begin
      npts=300 ; slow, but needed for decent convergence at small r
      dr=double(r[i])/npts
      rInterior = dindgen(npts)*dr
      sigmaInterior = nfw_sigma_offset(rInterior,roff[g],p,zl,r200=keyword_set(r200),r180=keyword_set(r180))
      integrand=rInterior*sigmaInterior ; sigmaInterior goes to infinity at r=0, so remove
      sigmaMean[i,g] = groupWeight[g] * 2./r[i]^2 * total(integrand[where(finite(integrand))]*dr)
   endfor
;stop
endfor

; \Delta\Sigma(R|R_{off})
deltaSigma = total((sigmaMean-sigmaR),2)

return, deltaSigma
end
