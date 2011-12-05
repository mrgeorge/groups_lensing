pro get_ds_model, fit_type, p_mean, zl, msun, x_mpc, $                             ; inputs
                  sis=sis, tis=tis, use_m200=use_m200,$                          ; switches
                  center=center,refcen=refcen,groupFile=groupFile,nfw_off=nfw_off,$; options for diff_stack
                  ps_term=ps_term, nfw_term=nfw_term, tot=tot,$                    ; output profiles
                  mnfw=mnfw,conc=conc,rnfw=rnfw                                    ; output NFW params

; determine model components for a given set of parameters
; this duplicates what ds_model.pro is supposed to do, but avoids the
; common block dependency for clarity

; fit_type is array telling which parameters are included in fit
; p_mean gives mean fit parameters
; zl and msun are the lens redshift and central stellar mass
; ps_term, nfw_term return model values at radii x_mpc
; if center, refcen, and groupFile are set, nfw_off will contain the
;   model NFW (assumed to be fit around refcen) convolved with the offset distribution between center
;   and refcen, to give the predicted signal around center
; note that nfw_term may include a centering offset if fit_type[6]>0

;-------------------------------------------------------------------------
; VARIABLES
;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
; Which parameters are fit [type = 1]  ?
; 0  M0     : baryonic mass
; 1  M_vir  : NFW virial MASS
; 2  C      : NFW concentration
; 3  alpha  : fraction
; 4  b      : bias
; 5  m_sigma: dispersion in central mass
; 6  offset : offset distance
;-------------------------------------------------------------------------

i=0

; M0 : 10^12 M_sun/h
if(fit_type[0] eq 1) then begin
   M0 = p_mean[i]
   i = i+1
endif else begin
; Fix from the data
   M0 = msun      ; alog10 units
endelse

; M_vir
if(fit_type[1] eq 1) then begin
   Mnfw = p_mean[i]             ; in LOG10
   i = i+1
endif 

; Concentration
if(fit_type[2] eq 1) then begin
   Conc=p_mean[i]
   i = i+1
endif else begin
   if(fit_type[2] EQ 2) then begin
      if (keyword_set(use_m200)) then begin
         if (keyword_set(use_maccio)) then begin
            conc = get_conc(Mnfw,zl,/use200,/maccio) ; Maccio
         endif else begin
            conc = get_conc(Mnfw,zl,/use200) ; Zhao
         endelse
      endif else begin
         conc = get_conc(Mnfw,zl,/virial)
      endelse
   endif
endelse

; q
if(fit_type[3] eq 1) then begin
   q = p_mean[i]
   i = i+1
endif else begin
   if(fit_type[3] eq 2) then begin ; Fix q
      ; Try a mass dependant model here :
      ; Base alpha on sm4 results

      alpha_sm4 = 0.4
      mcent_sm4 = 11.7
      
      alpha_0 = 0.0
      mcent_0 = 14.0            ; alpha is zero at these masses

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
endelse

; Bias
if(fit_type[4] eq 1) then begin
   bias = p_mean[i]
   i = i+1
endif else begin
   ; Fix from the data
   bias =  0.3 
endelse

   ; dispersion in mass
if(fit_type[5] eq 1) then begin
   m_sigma = p_mean[i]
   i = i+1
endif else begin
   ; Fix from the data
   m_sigma = 0.1 
endelse

; offset distance
if(fit_type[6] GT 0) then begin
   offset=p_mean[i]
   i=i+1
endif else begin
   offset=0.
endelse


; Calculate r200
if (keyword_set(use_m200) eq 1) then begin
   overdensity = get_overdensity( zl, /r200)
endif else begin
   overdensity = density_contrast(zl) ; virial
endelse
rho_crit = critical_density(zl)
factor   = alog10(double((4.0/3.0)*!Pi*overdensity*rho_crit))
r_log    = (1.0/3.0)*(Mnfw-factor)
rnfw     = 10.0^(r_log)


; Calculate model curves
if(n_elements(sis) GT 0) then begin                  ; replace point source with SIS
   r_eff=0.005 ; 5 kpc in Mpc
   meff_smeff_ratio=2. ; ratio of stellar mass within effective radius to total mass (DM+SM) within effective radius
   smeff_smtot_ratio=0.5 ; ratio of stellar mass within r_eff to "total" stellar mass 
   ; M0 ~ total stellar mass 
   ; M_eff = smtot * smeff/smtot * m_eff/smeff = DM+SM within r_eff
   m_eff=10.^(M0) * smeff_smtot_ratio * meff_smeff_ratio
   ps_term=m_eff/(2.*!pi*r_eff) / x_mpc / 1.e12 ; Msun/pc^2

endif else if(n_elements(tis) GT 0) then begin ; truncated isothermal PIEMD - see Kassiola & Kovner 1993 and Mira 2011
   r_eff=0.005 ; 5 kpc in Mpc
   r_cut=0.05 ; 50 kpc
   meff_smeff_ratio=2. ; ratio of stellar mass within effective radius to total mass (DM+SM) within effective radius
   smeff_smtot_ratio=0.5 ; ratio of stellar mass within r_eff to "total" stellar mass 
   ; M0 ~ total stellar mass 
   ; M_eff = smtot * smeff/smtot * m_eff/smeff = DM+SM within r_eff
   m_eff=10.^(M0) * smeff_smtot_ratio * meff_smeff_ratio ; Msun
   rho0=m_eff * (r_eff^2-r_cut^2) / (!pi * r_eff^2 * r_cut^2 * (!pi*r_eff - 4.*r_cut*atan(r_eff/r_cut))) ; Msun/Mpc^3

   ps_sigma_R=(rho0 * r_eff^2 * r_cut^2 * !pi)/(r_cut^2-r_eff^2) * (1./sqrt(r_eff^2 + x_mpc^2) - 1./sqrt(r_cut^2 + x_mpc^2)) ; Msun/Mpc^2 , Mira eq. 23.
   ps_sigmabar_R=(2.*rho0*r_eff^2*r_cut^2*!pi)/(x_mpc^2*(r_cut+r_eff)) * (1.-(sqrt(r_cut^2+x_mpc^2)-sqrt(r_eff^2+x_mpc^2))/(r_cut-r_eff)) ; Msun/Mpc^2 , Mira eq. 24.
   ps_term=(ps_sigmabar_R - ps_sigma_R) / 1.e12 ; Msun/pc^2
   
endif else $                             ; point source
   ps_term=10^(M0)/1.e12/(!pi*x_mpc^2); h^-1 Msun, factor of 1e12 to convert to pc^2

if(M0 EQ 0.) then ps_term[*]=0.

nfw_term=nfw_ds_offset(x_mpc,[rnfw,conc],zl,r200=keyword_set(use_m200),roff=offset)

tot=nfw_term + ps_term

if(keyword_set(center) AND keyword_set(refcen) AND keyword_set(groupFile)) then $
   nfw_off=nfw_ds_offset(x_mpc,[rnfw,conc],zl,groupFile,r200=keyword_set(use_m200),center=center,refcen=refcen)

end
