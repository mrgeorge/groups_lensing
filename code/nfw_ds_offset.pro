function get_proff, r, roff, offset_type

if(offset_type EQ 'delta3d') then begin
   proff=r/(roff^2 * sqrt(1.-r^2/roff^2))
endif else if(offset_type EQ 'max3d') then begin
   proff=(r/roff^2) * exp(-r^2/(2.*roff^2))
endif else begin
   print,'GET_PROFF: offset_type not recognized'
endelse

return, proff
end

function sigma_proff_integrand, x
; here x is a value for Roff, which is being integrated over
; roff_f is the parameter which describes the distribution P(Roff)

; return P(Roff) * \Sigma(R|Roff) for qromo integral to get \Sigma(R|P(Roff))
common common_integrand, p_f, zl_f, roff_f, mass_def_flag
common common_tabulate, reset_tab, reset_tab_proff, r_tab, nfw_sig_offset_tab, nfw_sig_proff_tab, this_r, offset_type

proff=get_proff(this_r,roff_f,offset_type)
if (mass_def_flag eq 0) then integrand=proff*tabulate_nfw_sigma_offset(this_r, x, p_f, zl_f,/r200)
if (mass_def_flag eq 1) then integrand=proff*tabulate_nfw_sigma_offset(this_r, x, p_f, zl_f,/r180)

return, integrand
end

function tabulate_nfw_sigma_proff,r, roff, p, zl, $
                                  r200=r200, r180=r180
; Similar to tabulate_nfw_sigma_offset but                                  
; Returns \Sigma(R|P(R_off))

common common_tabulate, reset_tab, reset_tab_proff, r_tab, nfw_sig_offset_tab, nfw_sig_proff_tab, this_r, offset_type

; Define variables if not previously defined
sz=size(r_tab)

if(sz[1] EQ 0 OR reset_tab_proff EQ 1) then begin ; Tabulate r in log space
   minRad=1.e-4 ; Mpc
   maxRad=4. ; Mpc
   npts=500 ; in log r
   step=(maxRad/minRad)^(1./(npts-1))
   r_tab=minRad*step^indgen(npts)
   nfw_sig_proff_tab = dindgen(npts) ; \Sigma(R|P(Roff))

   ; if offset_type='delta2d' then roff is simply a single projected 2d offset radius
   ; and \Sigma(R|P(Roff)) = \Sigma(R|Roff)
   ; other offset_type denote models for P(Roff) where roff is just a parameter
   ; and we need to integrate over P(Roff) to get \Sigma(R|P(R_off))
   if(offset_type EQ 'delta2d') then begin
      nfw_sig_proff_tab=tabulate_nfw_sigma_offset(r,roff,p,zl,r200=r200,r180=r180)
   endif else if(offset_type EQ 'delta3d') then begin ; roff is single 3d offset radius
      for ii=0,npts-1 do begin
         this_r=r_tab[ii]
         nfw_sig_proff_tab[ii]=qromo('sigma_proff_integrand',0.,roff,EPS=1.e-2)
      endfor
   endif else if(offset_type EQ 'max3d') then begin ; roff is 3d stddev in maxwellian
      for ii=0,npts-1 do begin
         this_r=r_tab[ii]
         nfw_sig_proff_tab[ii]=qromo('sigma_proff_integrand',0.,/midexp,EPS=1.e-2)
      endfor
   endif else begin
      print,'TABULATE_NFW_SIGMA_PROFF: unknown offset_type'
      stop
   endelse

   reset_tab_proff = 0
endif

; check limits on r
if(max(r) gt max(r_tab)) then begin
   print,'Stopped in tabulate_nfw_sigma_proff.pro, input r is out of limits'
   stop
endif

; Calculate using tabulated values
index=lonarr(n_elements(r))
res=fltarr(n_elements(r))
for ii=0,n_elements(r)-1 do index[ii]=min(where(r_tab GE r[ii]))

if(min(index) LT 0) then begin
   print,'Stopped in tabulate_nfw_sigma_proff.pro, invalid selection'
   stop
endif

; for r<minRad, fill with lowest tabulated ds value, else interpolate in log R
low=where(index EQ 0,nLow,complement=middle,ncomplement=nMiddle) ; low is where r<minRad, middle is where minRad<r<maxRad. r>maxRad has already been flagged
if(nLow GT 0) then res[low]=nfw_sig_proff_tab[0]
if(nMiddle GT 0) then res[middle]=interpol(nfw_sig_proff_tab,alog10(r_tab),alog10(r[middle]))

return, res
end

function tabulate_nfw_sigma_offset,r, roff, p, zl,$
                                  r200=r200, r180=r180

; Tabulate nfw_sigma_offset. Valid between 1e-4 and 4Mpc
; Could make smaller by using a smaller array
; Note: in this function input r can be a scalar or array

; Returns \Sigma(R|R_off)

common common_tabulate, reset_tab, reset_tab_proff, r_tab, nfw_sig_offset_tab, nfw_sig_proff_tab, this_r, offset_type


; Define variables if not previously defined
sz=size(r_tab)

if(sz[1] EQ 0 OR reset_tab EQ 1) then begin ; Tabulate r in log space
   minRad=1.e-4 ; Mpc
   maxRad=4. ; Mpc
   npts=500 ; in log r
   step=(maxRad/minRad)^(1./(npts-1))
   r_tab=minRad*step^indgen(npts)
   nfw_sig_offset_tab = dindgen(npts) ; \Sigma(R|Roff)

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

function nfw_rsigma_proff_integrand, x

common common_integrand, p_f, zl_f, roff_f, mass_def_flag

; return integrand = R*Sigma(R|P(Roff)) for determining mean Sigma(<R|P(Roff)) 
if (mass_def_flag eq 0) then res=x*tabulate_nfw_sigma_proff(x, roff_f, p_f, zl_f,/r200)
if (mass_def_flag eq 1) then res=x*tabulate_nfw_sigma_proff(x, roff_f, p_f, zl_f,/r180)

return, res
end

function nfw_ds_offset_sat, r, p, zl, roff, $
                            r200=r200, r180=r180, $
                            off_type=off_type

; To analytically calculate the satellite term in mock catalogs
; Alexie's version of Matt's code for gg-lensing + clustering project

; r is projected 2d distance in MPC
; p[0] is r200 the virial (or 200) radius in MPC
; p[1] is c the concentration parameter
; r180   -> 180 x mean density
; r200   -> 200 x critical density
; roff is projected 2d offset distance in Mpc

; Common block for the tabulation of nfw_sigma_offset
common common_tabulate, reset_tab, reset_tab_proff, r_tab, nfw_sig_offset_tab, nfw_sig_proff_tab, this_r, offset_type

reset_tab = 1
reset_tab_proff = 1
if(NOT(keyword_set(off_type))) then off_type='delta2d' ; unless specified, treat Roff as an explicit 2d radius
offset_type=off_type

; Common block for nfw_sigma_offset_integrand and nfw_sigma_proff_integrand (for qromo)
common common_integrand, p_f, zl_f, roff_f, mass_def_flag
p_f    = p
zl_f   = zl
roff_f = roff
if(keyword_set(r200)) then mass_def_flag=0
if(keyword_set(r180)) then mass_def_flag=1

; \Sigma(R|P(R_{off}))
sigmaR=dblarr(n_elements(r))
sigmaR=tabulate_nfw_sigma_proff(r,roff,p,zl,r200=r200,r180=r180)

; \bar{\Sigma}(<R|P(R_{off}))
sigmaMean=dblarr(n_elements(r))
sigmaMean=(2./r^2) * QROMO('nfw_rsigma_proff_integrand',replicate(1.e-4,n_elements(r)),r,EPS=1e-2) ; can increase precision here using EPS

; \Delta\Sigma(R|P(R_{off}))
ds=sigmaMean-sigmaR

return, ds
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
   deltaSigma+=groupWeight[g]*nfw_ds_offset_sat(r,p,zl,roff[g],off_type=off_type,r200=r200,r180=r180)
endfor

return, deltaSigma
end
