pro get_ds_model, fit_type, p_mean, lens_str, x_mpc, ps_term=ps_term, nfw_term=nfw_term,$
                  center=center,refcen=refcen,groupFile=groupFile,nfw_off=nfw_off,use_m200=use_m200,$
                  mnfw=mnfw,conc=conc,rnfw=rnfw
; determine model components for a given set of parameters
; this duplicates what ds_model.pro is supposed to do, but avoids the
; common block dependency for clarity

; fit_type is array telling which parameters are included in fit
; p_mean gives mean fit parameters
; lens_str is struct containing mean z and point source mass
; ps_term, nfw_term return model values at radii x_mpc
; if center, refcen, and groupFile are set, nfw_off will contain the
;   model NFW (assumed to be fit around refcen) convolved with the offset distribution between center
;   and refcen, to give the predicted signal around center
; note that nfw_term may include a centering offset if fit_type[6]>0

zl=lens_str.z_lens

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
   M0 = lens_str.msun_lens      ; alog10 units
endelse

; M_vir
if(fit_type[1] eq 1) then begin
   Mnfw = p_mean[i]             ; in LOG10
   i = i+1
endif 

; Concentration
if(fit_type[2] eq 1) then begin
   Conc = 10.0^(p_mean[i])           ; LOG
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
ps_term=10^(M0)/1.e12/(!pi*x_mpc^2) ; h^-1 Msun, factor of 1e12 to convert to pc^2
nfw_term=nfw_ds_offset(x_mpc,[rnfw,conc],zl,r200=keyword_set(use_m200),roff=offset)
if(keyword_set(center) AND keyword_set(refcen) AND keyword_set(groupFile)) then $
   nfw_off=nfw_ds_offset(x_mpc,[rnfw,conc],zl,groupFile,r200=keyword_set(use_m200),center=center,refcen=refcen)

end
