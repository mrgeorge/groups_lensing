function critical_density,z

;------------------------------------------------------------------------
;    
;   NAME         : critical_density : h^2 Msun Mpc-3
;
;   AUTHOR       : Alexie Leauthaud
;   DATE         : April  2007
;
;   PURPOSE      : Calculates a critical density at redshift z
;   RHO_CRIT = [3*h^2]/[8*Pi*Gn]
;-------------------------------------------------------------------------

defsysv,'!Omega_m', exists=exists
if NOT exists THEN define_cosmo

defsysv,'!Gn', exists=exists
if NOT exists THEN define_astro_cst

;--- h^2 Msun Mpc-3 ---
Hz    = hubble(z)             ; Hubble [h km s-1 Mpc-1]
rho_c = (Hz/!H0)^2*!Rho_c_0   ; h^2 Msun Mpc^-3
return,rho_c

end
