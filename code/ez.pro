function ez, z

; Hogg's useful E(z) function, Equation 14

defsysv,'!Omega_m', exists=exists
if NOT exists THEN define_cosmo

res = (!Omega_m*(1+z)^3+!Omega_k*(1+z)^2+!Omega_l)^(0.5)

return, res

end
