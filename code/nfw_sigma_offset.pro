function nfw_sigma_offset,r, roff, p, zl,$
                 r200 = r200, r180 = r180

; Following Johnston et al 2007 eq 32
; calls nfw_sigma and returns sigma(R|Roffset)
; see nfw_ds_offset to get delta_sigma(<R|Roffset)
; r is distance from given center in Mpc
; roff is projected 2d distance in Mpc between given center and true center of halo
; p[0] is virial radius in Mpc
; p[1] is c the concentration parameter
; zl is lens redshift

;npts=1000 ; doesn't matter much for small offsets, but important for precision when offsets are large
npts=max([30,30.*roff/0.01]) ; based on convergence tests, this should produce reliable profiles in the range of interest (~2% error at 50 kpc)
dtheta=2*!pi/npts
theta=dindgen(npts)*dtheta

arg=sqrt(abs(rebin([r^2],[n_elements(r),npts]) + roff^2 - 2.*roff*r#cos(theta))) ; offset radius, [nr x ntheta] array

; avoid arg=0 where r=roff
zero=where(arg EQ 0.,nZero)
if(nZero GT 0) then arg[zero]=min([arg[where(arg GT 0)],1.e-4]) ; assign arg=0 to a small positive value

sigma=reform(nfw_sigma(arg,p,zl,r200=r200,r180=r180),[n_elements(r),npts])

; assume random orientation of offsets and take azimuthal average
sigma_off=1./(2.*!pi) * total(sigma*dtheta,2)

return, sigma_off
end
