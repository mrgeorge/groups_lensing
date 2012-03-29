pro define_astro_cst

;--- CONSTANTS ---
long_Pi	=	3.141592653589
Pi_as	=	180.0*60*60		; Arcsec
Gn 	=	6.673e-11		; m3 kg-1 s-2
Log_Gn  =       -38.343                 ; log units of kpc3 s-2 M_solar-1 
C	= 	2.99792458e8		; m s-1
C_km	=	C/1000.0	 	; km s-1	

;--- CONVERSIONS ---
rad2deg	=	180.0/long_pi		; X_deg=X_rad*(180/Pi)
rad2as	=	(180.0/long_pi)*60*60	; X_as=X_rad*(180/Pi)*60*60

deg2rad	=	long_pi/180.0		;
as2rad	=	long_pi/(180.0*60*60)	; 

pc2m	=	3.08567758e16		; X_m=X_pc*pc2m
mpc2m	=	pc2m*1000000		; X_m=X_mpc*mpc2m
m2mpc	=	1.0/mpc2m

;--- SOLAR UNITS ---
Msun	=	1.989e33		; Solar mass in grammes

defsysv,'!Long_pi', exists=exists
IF NOT exists THEN defsysv,'!Long_pi',long_pi

defsysv,'!Gn', exists=exists
IF NOT exists THEN defsysv,'!Gn',gn

defsysv,'!mpc2m', exists=exists
IF NOT exists THEN defsysv,'!mpc2m',mpc2m

defsysv,'!as2rad', exists=exists
IF NOT exists THEN defsysv,'!as2rad',as2rad

defsysv,'!rad2as', exists=exists
IF NOT exists THEN defsysv,'!rad2as',rad2as


end
