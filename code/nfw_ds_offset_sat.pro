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
sigmaMean=(2./r^2) * QROMO('nfw_sigma_offset_function',replicate(1e-4,n_elements(r)),r,EPS=1e-2) ; can increase precision here using EPS

; \Delta\Sigma(R|R_{off})
res = sigmaMean-sigmaR

return, res
end
