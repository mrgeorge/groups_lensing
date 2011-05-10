pro plot_lensing_results, lensing_infile, ps_file, p_mean, fit_type,$
                          center=center,$
                          refcen=refcen,$
                          stackx=stackx,$
                          use_m200=use_m200,$
                          use_maccio=use_maccio, $
                          models=models

defsysv,'!Omega_m', exists=exists
if NOT exists THEN define_cosmo

;-------------------------------------------------------------------------
; Data Struct
;-------------------------------------------------------------------------
full_str=mrdfits(lensing_infile,1)
sel_str=where(full_str.e1_num GE 10)
zl = full_str.z_lens

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
;-------------------------------------------------------------------------

i=0

; M0 : 10^12 M_sun/h
if(fit_type[0] eq 1) then begin
   M0 = p_mean[i]
   i = i+1
endif else begin
; Fix from the data
   M0 = full_str.msun_lens      ; alog10 units
endelse

; M_vir
if(fit_type[1] eq 1) then begin
   Mnfw = p_mean[i]             ; in LOG10
   i = i+1
endif 

; Concentration
if(fit_type[2] eq 1) then begin
   Conc = 10.0^(P[i])           ; LOG
   i = i+1
endif else begin
   if(fit_type[2] EQ 2) then begin
      if (use_m200 eq 1) then begin
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

;-------------------------------------------------------------------------
; SET UP PLOT
;-------------------------------------------------------------------------
ps_open2,ps_file,/PORTRAIT,/ENCAPSULATED,XSIZE=7.0,YSIZE=5,/COLOR,THICK=3,/PS_FONT
device,/cmyk, /helvetica,font_size=14

; FILLED CIRCLE  
A = FINDGEN(17) * (!PI*2/16.)   
USERSYM, COS(A), SIN(A), /FILL
simpctable

xtitle=textoidl('Physical transverse distance,  R (h_{72}^{-1} Mpc)')
bar='!S!A=!R!N'
ytitle=textoidl('\Delta\Sigma = ')+bar+textoidl('\Sigma(<R) - \Sigma(R)  (h_{72} M')+sunsymbol()+textoidl(' pc^{-2})')

xst=1
yst=1
xlog=1
ylog=0
if(keyword_set(xlog)) then xr = [0.02,2] else xr=[0,1.5]
if(keyword_set(ylog)) then yr = [0.5,3000] else yr = [-200,400]
if(keyword_set(ylog)) then ytickf='loglabels' else ytickf=''
if(keyword_set(xlog)) then xtickf='loglabels' else xtickf=''

if(keyword_set(stackx)) then begin
   x    = full_str.plot_radius_kpc[sel_str]*rnfw
   y    = full_str.we1_mean[sel_str]
   yerr = full_str.we1_error[sel_str]
endif else begin
   x    = full_str.plot_radius_kpc[sel_str]/1e3
   y    = full_str.we1_mean[sel_str]
   yerr = full_str.we1_error[sel_str]
endelse

;-------------------------------------------------------------------------
; PLOT POINTS
;-------------------------------------------------------------------------
ploterror,x,y,yerr,xlog=xlog,ylog=ylog,yr=yr,xr=xr,xtickf=xtickf,ytickf=ytickf,$
          xtitle=xtitle,ytitle=ytitle,/xst,/yst,psym=8,symsize=1.0,color=!black


if(keyword_set(models)) then begin
   ;-------------------------------------------------------------------------
   ; PLOT THE MODEL
   ;-------------------------------------------------------------------------
   nxMpc = 50
   if(keyword_set(xlog)) then begin
      xbuffer=1.3
      x_mpc = 10.^(findgen(nxMpc)/(nxMpc-1)*alog10((xr[1]*xbuffer)/(xr[0]/xbuffer)))*xr[0]/xbuffer
   endif else x_mpc = findgen(nxMpc)/(nxMpc-1) * (xr[1]-xr[0]) + xr[0]
   
   ; Baryonic point source term
   if(fit_type[0] NE 0) then begin
      ps_term = 10^(M0)/1.e12/(!pi*x_mpc^2) ; h^-1 Msun, factor of 1e12 to convert to pc^2
      oplot,x_mpc,ps_term,color=!red,linestyle=1
   endif

   ; NFW term
   nfw = nfw_ds(x_mpc,[rnfw,conc],zl,r200=keyword_set(use_m200))
   oplot,x_mpc,nfw,color=!green,linestyle=2

   ; Offset NFW term
   if(keyword_set(center) AND keyword_set(refcen)) then begin
      nfw_off = nfw_ds_offset(x_mpc,[rnfw,conc],zl,r200=keyword_set(use_m200),center=center,refcen=refcen)
      oplot,x_mpc,nfw_off,color=!orange,linestyle=3
   endif

   ; Sum of terms (currently just NFW + point source if point source is included in model)
   if(fit_type[0] NE 0) then begin
      if(keyword_set(center) AND keyword_set(refcen)) then tot = ps_term + nfw_off $
      else tot = ps_term + nfw
      oplot,x_mpc,tot,color=!blue
   endif

   ;-------------------------------------------------------------------------
   ; Replot points
   ;-------------------------------------------------------------------------
   oploterror,x,y,yerr,psym=8,symsize=1.0,color=!black

   ;-------------------------------------------------------------------------
   ; Calculate chi^2
   ;-------------------------------------------------------------------------
   ; currently just NFW + point source if point source is included in model

   ; calculate expected model values at locations of data points
   if(keyword_set(center) AND keyword_set(refcen)) then $
      model_nfw = nfw_ds_offset(x,[rnfw,conc],zl,r200=keyword_set(use_m200),center=center,refcen=refcen) $
   else model_nfw = nfw_ds(x,[rnfw,conc],zl,r200=keyword_set(use_m200))

   if(fit_type[0] NE 0) then begin
      model_ps = 10^(M0)/1.e12/(!pi*x^2) ; h^-1 Msun, factor of 1e12 to convert to pc^2
      model_tot = model_ps + model_nfw
   endif else model_tot = model_nfw

   chisq = total((model_tot-y)^2/yerr^2)
   dof = n_elements(x)-n_elements(where(fit_type EQ 1))
endif

;-------------------------------------------------------------------------
; Legend
;-------------------------------------------------------------------------
nlens    = textoidl('N_{Lens}:')+string(full_str.lens,format="(I)")
z        = 'Redshift:'+string(full_str.z_lens,format="(f10.2)")
if(keyword_set(models)) then begin
   if keyword_set(use_m200) then begin
      m        = textoidl('log_{10}(M_{200}):')+string(mnfw,format="(f10.2)")
      r        = textoidl('R_{200}:')+string(rnfw,format="(f10.2)")
   endif else begin
      m        = textoidl('log(M_{vir}):')+string(mnfw,format="(f10.2)")
      r        = textoidl('R_{vir}:')+string(rnfw,format="(f10.2)")
   endelse
   c        = 'Concentration:'+string(Conc,format="(f10.2)")
   ;msun     = textoidl('log_{10}(M_{stellar}):')+string(stellar,format="(f10.2)") 
   ;eta      = 'SF eff:'+string(eta,format="(f10.1)") 
   ;bias     = 'Bias:'+string(bias,format="(f10.2)") 
   ;alpha    = 'Alpha:'+string(alpha,format="(f10.2)") 
   ;m_sigma  = 'm_sigma:'+string(m_sigma,format="(f10.1)") 

   ;box = flag to include/omit box around the legend (D=include)
   ;		  outline_color = color of box outline (D = !P.color)
   chisq_str = textoidl('\chi^2:')+string(chisq,format='(f10.2)')
   dof_str = 'd.o.f.:'+string(dof,format='(I)')
   items=[nlens,z,m,c,chisq_str,dof_str]
endif else items=[nlens,z]
legend,items,/right,linestyle=-99,box=0,spacing=1.5
    
ps_close

return
end
