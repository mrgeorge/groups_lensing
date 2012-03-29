function hubble, Z, use_a=use_a

;------------------------------------------------------------------------
;    
;   NAME         : hubble.pro
;
;   AUTHOR       : Alexie Leauthaud
;   DATE         : 31 Jan 2005
;
;   VERSION      : 1.0
;
;   PURPOSE      : Calculates Hubble constant at redshift z (km s-1 Mpc-1)
;
;   REQUIRES     : define comology constants: Omega_m, Omega_l
;
;   INPUT        : redshift
;   OUPUT        : H(Z) [100*h km s-1 Mpc-1]
;
;-------------------------------------------------------------------------

;-- This is for a flat universe (Omega_k = 0)
;-- H(z) = H0*E(Z) where E is introduced by D.Hogg

defsysv,'!Omega_m', exists=exists
if NOT exists THEN define_cosmo

; a = 1/(1+z)
if keyword_set(use_a) then begin
	h=!H0*sqrt(!Omega_m*(Z)^(-3)+!Omega_l)
endif else begin

    if (Z ge 0) then begin
	h=!H0*sqrt(!Omega_m*(1+Z)^(3)+!Omega_l)
    endif

endelse

return,h

end


