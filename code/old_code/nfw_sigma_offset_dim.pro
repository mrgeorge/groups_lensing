function nfw_sigma_offset_dim,r, roff, p, zl,$
                 r200 = r200, $
                 r180 = r180

; an attempt to avoid for loops with fitting NFW offsets using multi-dimensional arguments
; returns a 2d array of sigma values
;   x axis is r, y axis is roff (i.e. each group is a row, each radial bin is a column)

; Following Johnston et al 2007 eq 32
; calls nfw_sigma and returns sigma(R|Roffset)
; see nfw_ds_offset to get delta_sigma(<R|Roffset)
; r is distance from given center in Mpc
; roff is an ARRAY giving distance in Mpc between given centerS and true centerS of haloS
; p[0] is virial radius in Mpc
; p[1] is c the concentration parameter
; zl is lens redshift

npts=1000 ; doesn't matter much for small offsets, but important for precision when offsets are large
dtheta=2*!pi/npts
theta=dindgen(npts)*dtheta

; 3 vectors: r, roff, theta
; make one big matrix to hold all the Rprime values to send to nfw_sigma
nx=n_elements(r)
ny=n_elements(roff)
nz=npts
db3=systime(1)
m1=rebin(r^2,nx,ny,nz)
m2=rebin(reform(roff^2,1,ny),nx,ny,nz)
m3=2*rebin(r#roff,nx,ny,nz)*rebin(reform(cos(theta),1,1,nz),nx,ny,nz)
arg=sqrt(m1+m2+m3)
db4=systime(1)

arg_1d=reform(arg,nx*ny*nz)
db5=systime(1)
sigma_1d=nfw_sigma(arg_1d,p,zl,r200=keyword_set(r200),r180=keyword_set(r180))
db6=systime(1)
sigma=reform(sigma_1d,nx,ny,nz)
sigma_off=1./2*!pi * total(sigma*dtheta,3) ; average over third (theta) dimension, leaving sigma_off as a 2d array of r x roff sigma values

return, sigma_off
end
