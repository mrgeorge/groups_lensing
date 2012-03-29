pro define_cosmo

; Hubble Expansion Factor
little_h       =     0.72                  ; no unit 
H0	       =     100*little_h	   ; 100*h km s-1 Mpc-1

; Spaial Curvature : Flat Universe
Omega_k   =   0   

; Matter Energy Density (Omega_m = Omega_b + Omega_c + Omega_mu)
Omega_m   =  0.258

; Dark Energy
Omega_l	  =  1-(Omega_m+Omega_k)
wl        = -1

; Baryon Density (WMAP 5)
;wb        =  0.02273           ; Omega_b*h^2
;Omega_b   =  wb/(little_h^2)
Omega_b   =  0.0438

; Scalar spectral index at 0.0002/Mpc
ns        =     0.963

; Linear Theory Amplitude of matter. Fluctuations on 8 h^-1 Mpc
sigma_8   =     0.796

; Critical Density at z=0
; (This depends on h!)
Rho_c_0   =   critical_density_z0(little_h)   ; h^2 Msun Mpc^-3

; For h=1 -> 2.77536627*1e11        

;-------------------------------
print,'-----------------------------'
print,'Omega b  :', Omega_b
print,'Omega m  :', Omega_m
print,'Omega l  :', omega_l
;print,'little_h :', little_h
print,'spectral index',ns
print,'H0 :',H0
print,'Sigma_8  :', sigma_8
print,'Critical density, z=0 : ',Rho_c_0
print,'-----------------------------'
;-------------------------------

defsysv,'!H0', exists=exists
IF NOT exists THEN defsysv,'!H0',H0

defsysv,'!little_h', exists=exists
IF NOT exists THEN defsysv,'!little_h',little_h

defsysv,'!Omega_m', exists=exists
IF NOT exists THEN defsysv,'!Omega_m',Omega_m

defsysv,'!Omega_b', exists=exists
IF NOT exists THEN defsysv,'!Omega_b',Omega_b

defsysv,'!Omega_l', exists=exists
IF NOT exists THEN defsysv,'!Omega_l',Omega_l

defsysv,'!Omega_k', exists=exists
IF NOT exists THEN defsysv,'!Omega_k',Omega_k

defsysv,'!Rho_c_0', exists=exists
IF NOT exists THEN defsysv,'!Rho_c_0',Rho_c_0

defsysv,'!ns', exists=exists
IF NOT exists THEN defsysv,'!ns',ns

defsysv,'!sigma_8', exists=exists
IF NOT exists THEN defsysv,'!sigma_8',sigma_8

end
