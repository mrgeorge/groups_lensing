function get_omegam, z

;-- Get Omega_m at a given redshift --

;-- Redshift zero Cosmology --
defsysv,'!Omega_m', exists=exists
if NOT exists THEN define_cosmo

;-- Hogg EZ function --
Ez = ez(z)

;-- Evolution of Omega matter
res = (!Omega_m*(1+z)^3)/Ez^2   

RETURN,res

end
