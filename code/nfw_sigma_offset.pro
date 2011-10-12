function nfw_sigma_offset,r, roff, p, zl,$
                 r200 = r200, $
                 r180 = r180

; Following Johnston et al 2007 eq 32
; calls nfw_sigma and returns sigma(R|Roffset)
; see nfw_ds_offset to get delta_sigma(<R|Roffset)
; r is distance from given center in Mpc
; roff is distance in Mpc between given center and true center of halo
; p[0] is virial radius in Mpc
; p[1] is c the concentration parameter
; zl is lens redshift

npts=1000 ; doesn't matter much for small offsets, but important for precision when offsets are large
dtheta=2*!pi/npts
theta=dindgen(npts)*dtheta

arg=sqrt(rebin([r^2],[n_elements(r),npts]) + roff^2 - 2.*roff*r#cos(theta)) ; offset radius, [nr x ntheta] array
sigma=reform(nfw_sigma(arg,p,zl,r200=r200,r180=r180),[n_elements(r),npts])

; assume random orientation of offsets and take azimuthal average
sigma_off=1./(2.*!pi) * total(sigma*dtheta,2)

return, sigma_off
end
