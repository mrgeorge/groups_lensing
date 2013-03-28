function ds_integrand, rad
common sigma_profile, rkpc, sigma
thisSigma=interpol(sigma,alog10(rkpc),alog10(rad))
return, rad*thisSigma
end

function sigma_to_ds, rkpc, sigma

common sigma_profile, rad_kpc, sigma_msunkpc2
rad_kpc=rkpc
sigma_msunkpc2=sigma*1.e6 ; Msun/kpc^2 - assumes input sigma is in Msun/pc^2

sigmaMean=dblarr(n_elements(rkpc))
sigmaMean=(2./rkpc^2) * QROMO('ds_integrand',replicate(1.e-3,n_elements(rkpc)),rkpc,EPS=1e-2) ; can increase precision here using EPS

; \Delta\Sigma
ds=sigmaMean-sigma_msunkpc2 ; Msun / kpc^2

return, ds
end
