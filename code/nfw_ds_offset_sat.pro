function nfw_ds_offset_sat, r, p, zl, roff, $
                        r200 = r200, $
                        r180 = r180, $
                        sigma=sigma

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
for i=0,n_elements(r)-1 do begin
   sigmaR[i] = tabulate_nfw_sigma_offset(r[i],roff,p,zl,r200=keyword_set(r200),r180=keyword_set(r180))
endfor

; Just calculate SIGMA :
if keyword_set(sigma) then begin
   res = sigmaR

; Calculate DeltaSigma :
endif else begin

   ; \bar{\Sigma}(<R|R_{off})
   for i=0,n_elements(r)-1 do begin
      sigmaMean[i] = (2./r[i]^2) * QROMO('nfw_sigma_offset_function', 1e-4, r[i] , /DOUBLE , EPS=1e-2 ) ; can increase precision here using EPS
   endfor

   ; \Delta\Sigma(R|R_{off})
   res = sigmaMean-sigmaR
endelse

return, res
end
