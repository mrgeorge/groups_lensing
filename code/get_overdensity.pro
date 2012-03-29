function get_overdensity, zl, r200=r200, r180=r180, cr180=cr180,r500=r500, b200=b200

;-------------------------------------------------------------------------
; PURPOSE:
;  Get the overdensity defined wrspt to the CRITICAL density
;  Default is virial overdensity
;-------------------------------------------------------------------------

;-- Define Cosmology --
defsysv,'!Omega_m', exists=exists
if NOT exists THEN define_cosmo

;-- Virial overdensity --
overdensity = density_contrast(zl)

;-- Overdensity of 200 --
if keyword_set(r200) then overdensity = 200.0

; -- Overdensity of 180c
if keyword_set(cr180) then overdensity = 180.0

; -- 180 * rho_mean : define wspt the critical density here --
; 180 x background = 180 x rho_mean(z) = 180 x Omega_m(z) x rho_crit(z)
if keyword_set(r180) then overdensity = 180 * get_omegam(zl)

; -- 200 times background --
if keyword_set(b200) then overdensity = 200 * get_omegam(zl)

;-- Overdensity of 200 --
if keyword_set(r500) then overdensity = 500.0

; -- Fixed radius : need to divide by rho_crit
;if keyword_set(r500fix) then overdensity = 500/critical_density(zl)

RETURN,overdensity

end
