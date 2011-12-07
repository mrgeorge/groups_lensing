function nfw_ds_offset_dim, r, p, zl, $
                        r200 = r200, $
                        r180 = r180, $
                        center=center, $
                        refcen=refcen 

; an attempt to avoid for loops with fitting NFW offsets using multi-dimensional arguments

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

; get array off offset lengths
groupsel=group[where(group.flag_include EQ 1)]
get_center_coords,groupsel,center,ra1,dec1,sel1
get_center_coords,groupsel,refcen,ra2,dec2,sel2
roff=distance(ra1,dec1,ra2,dec2)*3600.*groupsel.lensing_r200_mpc/groupsel.lensing_r200_as ; Mpc
offset_kpc=roff*1000.
max_offset_kpc=50.
sel=where(ra1 GT 0 AND ra2 GT 0 AND offset_kpc LT max_offset_kpc)
roff=roff[sel]
ngroups=n_elements(sel)

; get relative weights from number of sources in bin behind lens
; (need to modify gglensing_v6 to get this - use 1 for now)
groupWeight=replicate(1./ngroups,ngroups)

; \Sigma(R|R_{off})
; n_elements(r) x ngroups matrix with sigma(r) values for each group
sigmaR = nfw_sigma_offset_dim(r,roff,p,zl,r200=keyword_set(r200),r180=keyword_set(r180))

; \bar{\Sigma}(<R|R_{off})
sigmaMean=dblarr(n_elements(r),ngroups)
npts=300 ; slow, but needed for decent convergence at small r
for i=0, n_elements(r)-1 do begin
   dr=double(r[i])/npts
   rInterior= dindgen(npts)*dr
   sigmaInterior = nfw_sigma_offset_dim(rInterior,roff,p,zl,r200=keyword_set(r200),r180=keyword_set(r180))
   integrand = rebin(rInterior,npts,ngroups)*sigmaInterior
   sigmaMean[i,*] = 2./r[i]^2 * total(integrand[where(finite(integrand))]*dr,1)
endfor

; \Delta\Sigma(R|R_{off})
; collapse matrix to n_elements(r) array, with group weightings
deltaSigma = total((sigmaMean-sigmaR)*rebin(reform(groupWeight,1,ngroups),n_elements(r),ngroups),2)

return, deltaSigma
end
