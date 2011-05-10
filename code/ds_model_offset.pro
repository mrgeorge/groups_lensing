function ds_model_offset, X, P,$
                   nfw=nfw,$
                   point_s=point_s,$
                   group=group,$
                   two_halo=two_halo,$
                   weak_nfw=weak_nfw,$
                   ds_z0=ds_z0,$
                   return_all=return_all,$
                   skip_dslin=skip_dslin

; NOV 2008 -> CHANGED TO USE MASS AS A VARIABLE INSTEAD OF R200 !!!
; Get_conc now takes m in alog10 units

;------------------------------------------------------------------------
;   FUNCTION TO DESCRIBE DELTA SIGMA
;
;   INPUT        : X	Array radius to calculate shear !! in Mpc !!
;		   P	Array of model values 
;		   
;   OUPUT        : array of delta sigma
;
;   NOTE         : Overdensity is Virial
;-------------------------------------------------------------------------

defsysv,'!Omega_m', exists=exists
if NOT exists THEN define_cosmo

;start=systime(1)

;-------------------------------------------------------------------------
; Common Block
;-------------------------------------------------------------------------

common fit_options, q_c, lens_redshift, fit_type, lens_m_sun, log_sm, use_m200, neg_points, pos_points,str2,str3, use_group, use_maccio, xmar, ymar, xchars, ychars, no_title, ws_corr, lz_mean, sx
common xi_lin_variable, dslin   ; Just for this program
common xGroupCenter, xCenter, xRefCen

zl = lens_redshift

;-------------------------------------------------------------------------
; Do we ignore certain componants [fit = 0]? 
; 0  baryonic mass
; 1  M_vir : NFW MASS
; 2  C     : NFW concentration
; 3  alpha : fraction
; 4  Bias
; 5  m_sigma : dispersion in central mass
;-------------------------------------------------------------------------

b1 = 1.0 & n1 = 1.0 & g1 = 1.0 & t1 = 1.0 & s1 = 1.0

; Baryonic term
if (fit_type[0] eq 0) then begin
    b1 = 0.0
endif

; NFW
if ((fit_type[1] eq 0) or (fit_type[2] eq 0)) then begin 
    n1 = 0.0
endif

; Group
if ((fit_type[3] eq 0)) then begin
    alpha = 0.0
    g1 = 0.0
endif

; Two halo
if (fit_type[4] eq 0) then begin 
    t1 = 0.0
    bias = 0.0
endif

; Dispersion in central mass
if (fit_type[5] eq 0) then begin 
    s1 = 0.0 
endif

; If we ignore a term then return array of zero's
zero_array    = fltarr(n_elements(X))
zero_array(*) = 0.0

;-------------------------------------------------------------------------
; Which parameters are fit [type = 1]  ?
; 0  M0     : baryonic mass
; 1  M_vir  : NFW virial MASS
; 2  C      : NFW concentration
; 3  alpha  : fraction
; 4  b      : bias
; 5  m_sigma: dispersion in central mass
;-------------------------------------------------------------------------

i = 0

; M0 : 10^12 M_sun/h
if(fit_type[0] eq 1) then begin
    M0 = P[i]
    i = i+1
endif else begin
    ; Fix from the data
    M0 = log_sm   ; alog10 units
endelse

; R_vir : Mpc/h
if(fit_type[1] eq 1) then begin
    Mnfw = P[i]  ; in LOG10
    i = i+1
endif 

; Concentration
if(fit_type[2] eq 1) then begin
    Conc = 10.0^(P[i])      ; LOG
    i = i+1
endif 

; q
if(fit_type[3] eq 1) then begin
    q = P[i]
    i = i+1
endif
if(fit_type[3] eq 2) then begin ; Fix q
    ; Try a mass dependant model here :
    ; Base alpha on sm4 results

    alpha_sm4 = 0.4
    mcent_sm4 = 11.7

    alpha_0 = 0.0
    mcent_0 = 14.0              ; alpha is zero at these masses

    alpha_a = (alpha_sm4-alpha_0)/(Mnfw_sm4-mcent_0) ; linear function in log space
    alpha_b = alpha_sm4-(alpha_a*mcent_sm4)

    alpha = (alpha_a*Mnfw)+alpha_b

    ; don't let alpha get too big
    if(alpha gt 0.4) then alpha=0.4

    ; make sure doesn't go neg or zero
    if (Mnfw gt 13.9) then alpha=0.01

    beta   =  alpha/(1-alpha) 
    q      =  alog(beta)

endif 

; Bias
if(fit_type[4] eq 1) then begin
    bias = P[i]
    i = i+1
endif else begin
    ; Fix from the data
    bias =  0.3 
endelse

; dispersion in mass
if(fit_type[5] eq 1) then begin
    m_sigma = P[i]
    i = i+1
endif else begin
    ; Fix from the data
    m_sigma = 0.1 
endelse

;-------------------------------------------------------------------------
; ********************* NFW TERM *****************************************
;-------------------------------------------------------------------------

if (n1 ne 0) then begin

; Is Concentration is fixed by m-c relation ?
if(fit_type[2] eq 2) then begin
    if (use_m200 eq 1) then begin
        if (use_maccio eq 1) then begin
            Conc = get_conc(Mnfw,zl,/use200,/maccio)  ; Maccio
        endif else begin
            Conc = get_conc(Mnfw,zl,/use200)          ; Zhao
        endelse
    endif else begin
        Conc = get_conc(Mnfw,zl,/virial)
    endelse
endif

; Calculate r200
if (use_m200 eq 1) then begin
    overdensity = get_overdensity( zl, /r200)
endif else begin
    overdensity = density_contrast(zl) ; virial
endelse

rho_crit = critical_density(zl)
factor   = alog10(double((4.0/3.0)*!Pi*overdensity*rho_crit))
r_log    = (1.0/3.0)*(Mnfw-factor)
rnfw     = 10.0^(r_log)

if (sx gt 0) then begin
    X2=X*rnfw*1e3                ;1e3 to compensate for the division by 1e3 earlier
endif else begin
    X2=X
endelse

; Include the dispersion in mass selection here
if (fit_type[5] eq 0) then begin 
    if (use_m200 eq 1) then begin
        nfw_term = nfw_ds(X2,[rnfw,Conc],zl,/r200)
;        nfw_term = nfw_ds_offset_dim(X2,[rnfw,Conc],zl,/r200,center=xCenter,refcen=xRefCen)
    endif else begin
        nfw_term = nfw_ds(X2,[rnfw,Conc],zl)
;        nfw_term = nfw_ds_offset_dim(X2,[rnfw,Conc],zl,center=xCenter,refcen=xRefCen)
    endelse
endif else begin
; This is with a central dispersion
; mnfw is not in log units
    nfw_term = nfw_ds_sig( X2, 10.0^(Mnfw), Conc, m_sigma, zl)
print,'>>>'
stop
endelse


; nfw term for non weak shear
if (ws_corr eq 1) then begin

    ; get sigma
       if (use_m200 eq 1) then begin
        sig = nfw_sig(X2,[rnfw,Conc],zl,/r200)
    endif else begin
        sig = nfw_sig(X2,[rnfw,Conc],zl)
    endelse

    ; Calculate the ws term: Sigma * DeltaSigma * Lz
    nfw_ws_term = sig*nfw_term*lz_mean

endif else begin
    nfw_ws_term = zero_array
endelse


endif else begin
    nfw_term = zero_array
endelse

;-------------------------------------------------------------------------
; ********************* POINT SOURCE TERM ********************************
;-------------------------------------------------------------------------

if (b1 ne 0) then begin

;-----------------------------------------
; M0    [ h^-1 Msun  ]
; X^2   [ h^-2 Mpc^2 ]   -> h^1 Msun pc^-2
;-----------------------------------------

; Convert Point Source
M0 = 10.0^(M0)   ; h^-1 Msun

ps_term  = ( M0/1e12 )/(!pi*X2^2)      ; factor of 1e12 to convert to pc^2

endif else begin
    ps_term = zero_array
endelse

;-------------------------------------------------------------------------
; ********************* NFW TERM AT Z=0 **********************************
;-------------------------------------------------------------------------

if (n1 ne 0) then begin

    ; HOW TO DO THIS SCALING ??

; convert to rs and delta

    ;if (use_m200 eq 1) then begin
    ;    res = convert_nfw(mnfw, 0.0, 0.0, Conc, 0.0, lens_redshift, /r200)
    ;    rs    = res[2] 
    ;    delta = res[4]
    ;    print,''
    ;    print,'mass',mnfw
    ;    print,'rnfw',rnfw
    ;    print,'rs',rs
    ;    print,'delta',delta
    ;    print,'conc',Conc
    ;endif else begin
    
    ;endelse

; convert back to rnfw and Conc at z=0
    ;if (use_m200 eq 1) then begin
    ;    res = convert_nfw(0.0, 0.0, rs, 0.0, delta, 0.0, /r200)
    ;    mnfw_0    = res[0] 
    ;    rnfw_0    = res[1]
    ;    conc_0    = res[3]
    ;    print,''
    ;    print,'mass 0',mnfw_0
    ;    print,'rnfw_0',rnfw_0
    ;    print,'conc_0',conc_0
    ;endif else begin
    
    ;endelse

    ; == Calculate critical density
    ;rho_1 = critical_density(lens_redshift)
    ;rho_0 = critical_density(0.0)

    ;if (use_m200 eq 1) then begin
    ;    nfw_term_z0 = nfw_ds(X2,[rnfw,Conc],zl,/r200)*(rho_0/rho_1)
    ;endif else begin
    ;    nfw_term_z0 = nfw_ds(X2,[rnfw,Conc],zl)*(rho_0/rho_1)
    ;endelse

    ; How would the halo of the same mass look like at z=0?
    if (use_m200 eq 1) then begin
        nfw_term_z0 = nfw_ds(X2,[rnfw,4.0],0.0,/r200)
    endif else begin
        nfw_term_z0 = nfw_ds(X2,[rnfw,4.0],0.0)
    endelse

endif else begin
    nfw_term_z0 = zero_array
endelse

;-------------------------------------------------------------------------
; ********************* GROUP TERM ***************************************
;-------------------------------------------------------------------------

if (g1 ne 0) then begin

   ; ADD THE MASS HERE !!!!!
   bump_term = get_ds_1halo_sat_interpol(X2,alog10(mnfw))

   ; Fraction
   alpha = 1/(1+exp(-q))  

   ; Group term
   group_term = alpha*bump_term

endif else begin
    group_term = zero_array
endelse

;-------------------------------------------------------------------------
; ********************* TWO HALO TERM ************************************
;-------------------------------------------------------------------------

if (t1 ne 0) then begin

    ; This part only needs to be calculated once (called once in ds_model.pro)
    ; NOTE ; There is a floating underflow problem here that is unsolved...
    ; (it comes from linear_xi, probably not so important...)

    if NOT keyword_set(skip_dslin) then begin
          ; Calculate linear_ds (assumes PHYSICAL distances)
          lin_ds, X2, zl, dslin
    endif

    ;-------------------------------------------------------------
    ; Bias
    ; NOTE : The mass must be in M180b in order to derive the bias
    ;-------------------------------------------------------------
    
    ; CONVERSION PROBLEM HERE ??? hummm...I don't remember what the problem is...
    if (fit_type[4] eq 2) then begin

        if (use_m200 eq 1) then begin
            mm = convert_nfw_masses(10.0^mnfw,Conc,zl,/m200_2_m180)    
        endif else begin
            mm = convert_nfw_masses(10.0^mnfw,Conc,zl,/mvir_2_m180)
        endelse
        
        m180 = mm[0]    
 
        ; Mass must be in log units
        bias = mass_2_bias(alog10(m180),zl)

        if(alog10(m180) lt 10 or alog10(m180) gt 17) then begin
            print,'stopped in ds_model.pro'
            print,'m180 has funny value when calculating the bias'
            stop
        endif
    endif

    ;---------------------------------------------------------------
    ; Two halo term in physical units
    ;---------------------------------------------------------------

    two_halo_term = bias*dslin

endif else begin
    two_halo_term = zero_array
endelse

;-------------------------------------------------------------------------
; ********************* TOTAL MODEL ************************************
;-------------------------------------------------------------------------

ds = ps_term + nfw_term + nfw_ws_term+ group_term + two_halo_term

;-------------------------------------------------------------------------
; Return only a given term (to make plots)
;-------------------------------------------------------------------------

; Return this small number if we are supposed to ignore this term
fake_it = fltarr(n_elements(X2))
fake_it(*) = 1e-3

if keyword_set(nfw) then begin
    ds = nfw_term
    if (n1 eq 0) then ds=fake_it
endif
 
if keyword_set(weak_nfw) then begin
    ds = nfw_ws_term
endif

if keyword_set(ds_z0) then begin
    ds = nfw_term_z0
    if (n1 eq 0) then ds=fake_it
endif

if keyword_set(point_s) then begin
    ds = ps_term
    if (b1 eq 0) then ds=fake_it
endif

if keyword_set(group) then begin
    ds = group_term
    if (g1 eq 0) then ds=fake_it
endif

if keyword_set(two_halo) then begin
    ds = two_halo_term
    if (t1 eq 0) then ds=fake_it
endif

if keyword_set(return_all) then begin
    PP = [ M0, Rnfw, Conc, alpha, bias, mnfw, m_sigma ]
    ds = PP
endif

;print,''
;print,'Runtime = ',(systime(1)-start),' sec' 
;stop

;-------------------------------------------------------------------------
; THE END
;-------------------------------------------------------------------------

RETURN,ds

end
