function rho_interp, r3d, rho, thisR
; interpolate 3d density profile to thisR
thisRho=interpol(rho,alog10(r3d),alog10(thisR)) ; assuming r3d are log-spaced, do linear interpolation in log(rad)
return, thisRho>0 ; restrict density >=0
end

function sigma_integrand, theta
; projection from rho to sigma is an integral over theta
; use common block to get rho profile at r3d=R2d/cos(theta)

common density_profile, r3d, rho, thisR2d

thisRho=rho_interp(r3d,rho,thisR2d/cos(theta))
integrand=thisRho * thisR2d/cos(theta)^2

return, integrand
end

function rho_to_sigma, rho3d, rkpc3d, rkpc2d
; take a tabulated 3d density profile rho(r3d)
; interpolate as needed to sample at arbitrary points
; integrate along LOS to get projected density profile sigma(r2d)

common density_profile, r3d, rho, thisR2d
r3d=rkpc3d
rho=rho3d

nout=n_elements(rkpc2d)
sigma=fltarr(nout)

for ii=0,nout-1 do begin
   thisR2d=rkpc2d[ii] ; set common block variable
   sigma[ii]=2.*qromo('sigma_integrand',0.,!pi/2)
endfor

return,sigma
end
