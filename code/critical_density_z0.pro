function critical_density_z0, h

;------------------------------------------------------------------------
;    
;   NAME         : critical_density : h^2 Msun Mpc-3
;
;   AUTHOR       : Alexie Leauthaud
;   DATE         : April  2007
;
;   PURPOSE      : Calculates a critical density at z=0
;   RHO_CRIT = [3*h^2]/[8*Pi*Gn]
;-------------------------------------------------------------------------

x = 277.536627*(h)^2    ; Msun Kpc^-3
x = x * 1e9             ; Msun Mpc^-3

return,x

end
