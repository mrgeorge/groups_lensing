function  get_conc, m, z, virial=virial, use200=use200, maccio=maccio, munoz=munoz

;---------------------------------------------------------------------------------
;  WMAP 5 cosmology : O_m = 0.258, h=0.72, sigma_8=0.796, ns=0.963,Omega_b=0.0438
;  Default -> Zhao et al (only implemented Zhao)
;  Opional -> Maccio et al (200 and virial)
;             (For ALL halos (not just relaxed): Fig 2 and 3 in paper)
;  Uses M200 and C200
;  M200 IS NOW IN LOG UNITS !!!!
;  NOTE : M200 input uses h=0.72
;  NOTE: munoz comes from MUNOZ et al. 2010 and the text file has been calculated in
;  mcmc_lensing.c and mass is in units of 200 * background
;---------------------------------------------------------------------------------

defsysv,'!Omega_m', exists=exists
if NOT exists THEN define_cosmo

if (m gt 50) then begin
    print,'stopped in get_conc.pro -> now uses log10 mass as an input'
    stop
endif

;-------------------------------------------------------------------------
; >>>  ZHAO
;-------------------------------------------------------------------------

if (NOT keyword_set(maccio) and NOT keyword_set(munoz)) then begin

    c_z = get_zhao_conc(m, z)   ; this has already accounted for h & z

    if keyword_set(virial) then begin
        print,'I have not yet calcualte virial for Zhao'
        print,'stopping in get_conc.pro'
        stop
    endif
endif

;-------------------------------------------------------------------------
; >>>  MACCIO
;      According to Jeremy: h=1 in these masses
;-------------------------------------------------------------------------

if (keyword_set(maccio)) then begin

    if keyword_set(use200) then begin
        ; EQ 10 Maccio et al.
        ; mass_0.72 = mass_1/!little_h 
        ; mass_0.72 * !little_h =  mass_1 
        log_c_0 = 0.787-0.110*(m + alog10(!little_h)-12)
        c_0=10.0^(log_c_0)
        c_z = c_0 *(ez(z))^(-2.0/3.0) ; eq 19 Maccio et al.
    endif

    if keyword_set(virial) then begin
        ; EQ 9 Maccio et al.
        log_c_0 = 0.925-0.108*(m + alog10(!little_h)-12)
        c_0=10.0^(log_c_0)
        ;print,'in get_conc.pro need to define he z evolution here'
        ;print,'see equation 17 in maccio et al'
        ;stop
        ;c_z = c_0/(1+z)
        c_z = c_0/(1+z)
    endif
endif

;-------------------------------------------------------------------------
; >>>  MUNOZ
; back200 masses
; masses assumed in get_munoz_conc are in units of h^-1
;-------------------------------------------------------------------------

if (keyword_set(munoz)) then begin
   ; need to transform into h masses
   c_z = get_munoz_conc(m + alog10(!little_h), z) 
endif

;-------------------------------------------------------------------------
;-------------------------------------------------------------------------

if ((not keyword_set(use200)) and (not keyword_set(virial)) and (not keyword_set(munoz))) then begin
    print,'please specify which concentration to use in get_conc.pro'
    print,'stopping in get_conc.pro'
    stop
endif

if (c_z le 0.1 or c_z gt 25) then begin
    print,'stopped in get_conc.pro -> c_z has an unusually high or low value ...should check'
    print,'CONC: ',c_z
    stop
endif

; Redshift dependance:
return,c_z
end
