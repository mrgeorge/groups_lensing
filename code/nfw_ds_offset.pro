function nfw_sigma_offset_integrand, x

common common_nfw_sigma_offset_function, p_f, zl_f, roff_f, mass_def_flag

; return integrand = R*Sigma(R|Roff) for determining mean sigma(<R|roff) 
if (mass_def_flag eq 0) then res=x*tabulate_nfw_sigma_offset(x, roff_f, p_f, zl_f,/r200)
if (mass_def_flag eq 1) then res=x*tabulate_nfw_sigma_offset(x, roff_f, p_f, zl_f,/r180)

return, res
end

function tabulate_nfw_sigma_offset,r, roff, p, zl,$
                        r200 = r200, $
                        r180 = r180

; Tabulate nfw_sigma_offset. Valid between 1e-4 and 4Mpc
; Could make smaller by using a smaller array
; Note: in this function input r can be a scalar or array

common common_nfw_sig_offset, reset_tab, r_tab, nfw_sig_offset_tab 

; Define variables if not previously defined
sz=size(r_tab)

if(sz[1] EQ 0 OR reset_tab EQ 1) then begin ; Tabulate r in log space
   minRad=1.e-4 ; Mpc
   maxRad=4. ; Mpc
   npts=500 ; in log r
   step=(maxRad/minRad)^(1./(npts-1))
   r_tab=minRad*step^indgen(npts)

   nfw_sig_offset_tab = dindgen(npts)
   nfw_sig_offset_tab=nfw_sigma_offset(r_tab,roff,p,zl,r200=r200,r180=r180)
   reset_tab = 0
endif

; check limits on r
if(max(r) gt max(r_tab)) then begin
   print,'Stopped in tabulate_nfw_sigma_offset.pro, input r is out of limits'
   stop
endif

; Calculate using tabulated values
index=lonarr(n_elements(r))
res=fltarr(n_elements(r))
for ii=0,n_elements(r)-1 do index[ii]=min(where(r_tab GE r[ii]))

if(min(index) LT 0) then begin
   print,'Stopped in tabulate_nfw_sigma_offset.pro, invalid selection'
   stop
endif

; for r<minRad, fill with lowest tabulated ds value, else interpolate in log R
low=where(index EQ 0,nLow,complement=middle,ncomplement=nMiddle) ; low is where r<minRad, middle is where minRad<r<maxRad. r>maxRad has already been flagged
if(nLow GT 0) then res[low]=nfw_sig_offset_tab[0]
if(nMiddle GT 0) then res[middle]=interpol(nfw_sig_offset_tab,alog10(r_tab),alog10(r[middle]))

return, res
end

function nfw_ds_offset_sat, r, p, zl, roff, $
                        r200 = r200, $
                        r180 = r180

; To analytically calculate the satellite term in mock catalogs
; Alexie's version of Matt's code for gg-lensing + clustering project

; r is distance in MPC
; p[0] is r200 the virial (or 200) radius in MPC
; p[1] is c the concentration parameter
; r180   -> 180 x mean density
; r200   -> 200 x critical density
; roff is offset distance in Mpc

; Common block for the tabulation of nfw_sigma_offset
common common_nfw_sig_offset, reset_tab, r_tab, nfw_sig_offset_tab 
reset_tab = 1

; Common block for nfw_sig_function (for qromo)
common common_nfw_sigma_offset_function, p_f, zl_f, roff_f, mass_def_flag

p_f    = p
zl_f   = zl
roff_f = roff

if(keyword_set(r200)) then mass_def_flag=0
if(keyword_set(r180)) then mass_def_flag=1

; \Sigma(R|R_{off})
sigmaR    = dblarr(n_elements(r))
sigmaMean = dblarr(n_elements(r))

; r is an array :
sigmaR=tabulate_nfw_sigma_offset(r,roff,p,zl,r200=keyword_set(r200),r180=keyword_set(r180))

; \bar{\Sigma}(<R|R_{off})
sigmaMean=(2./r^2) * QROMO('nfw_sigma_offset_integrand',replicate(1e-4,n_elements(r)),r,EPS=1e-2) ; can increase precision here using EPS

; \Delta\Sigma(R|R_{off})
res = sigmaMean-sigmaR

return, res
end


function nfw_ds_offset, r, p, zl, groupFile, $
                        r200 = r200, r180 = r180, $
                        center=center, refcen=refcen, $
                        roff=roff, off_type=off_type

; center is the string name of the center under consideration
; refcen is the best halo center, default 'mmgg_scale'
; roff can be used as a free parameter for offset fitting
;    its behavior is determined by off_type:
;       'delta2d'-> offset distribution is a delta function in 2d projected radius
;                   roff is that 2d radius
;       'delta3d'-> offset distribution is a delta function in 3d radius
;                   roff is that 3d radius
;       'max3d'-> offset distribution is a Maxwellian in 3d radius
;                   roff is the dispersion parameter (=std.dev. of offset in each dimension)
; roff and off_type set the offset (distribution) explicitly. If not
; set, the offset distribution is taken to be the distribution between center and refcen
; r is distance in MPC to be modeled
; p[0] is r200 the virial (or 200) radius in MPC
; p[1] is c the concentration parameter
; r180   -> 180 x mean density
; r200   -> 200 x critical density

if(n_elements(roff) GT 0) then begin ; SPECIFY OFFSETS DIRECTLY
   if(min(roff) LT 0.) then begin
      print,'nfw_ds_offset returning ds=0: roff=',roff
      return, 0.*r              ; return 0's since roff shouldn't be negative
   endif
   if(max(roff) EQ 0.) then return, nfw_ds(r,p,zl,r200=r200,r180=r180) ; no offsets

   ; (else)
   ngroups=n_elements(roff)
endif else begin ; USE OFFSETS BETWEEN CENTERS FROM GROUP CATALOG
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
