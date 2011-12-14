function nfw_rho,r,p,zl
;p[0] is r200 the virial radius in Mpc
;p[1] is c the concentration parameter

; copied from /Users/alexie/idl/dave_idl/funcs/
; with updates using Alexie's code
;  assumes overdensity=200.


defsysv,'!Omega_m', exists=exists
if NOT exists THEN define_cosmo

r200=p[0]
conc=p[1]
x=r/r200

rho_crit = critical_density(zl) ; Msun Mpc^-3

del_c=(200./3.0)*conc^3 /(alog(1.0+conc)-conc/(1.0+conc))

return, rho_crit*del_c /(conc*x * (1.0+conc*x)^2)
end
