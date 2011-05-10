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

; loop over each element of r
sigma_off=dblarr(n_elements(r))
for i=0,n_elements(r)-1 do begin

   ;offset radius - array of length npts(theta)
   arg=sqrt(r[i]^2 + roff^2 + 2.*r[i]*roff*cos(theta))

   sigma=nfw_sigma(arg,p,zl,r200=keyword_set(r200),r180=keyword_set(r180))

   ; assume random orientation of offsets and take azimuthal average
   sigma_off[i]=1./(2*!pi) * total(sigma*dtheta)
endfor

return, sigma_off
end
