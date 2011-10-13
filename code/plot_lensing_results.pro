pro plot_lensing_results, lensing_infile, ps_file, p_mean, fit_type,$
                          center=center,$
                          refcen=refcen,$
                          groupFile=groupFile,$
                          stackx=stackx,$
                          use_m200=use_m200,$
                          use_maccio=use_maccio, $
                          models=models,$
                          fit_type2=fit_type2,$
                          p_mean2=p_mean2

; plot delta sigma vs r from lensing infile
; if models keyword is set, overplot model curves using fitted parameters p_mean
;   based on what parameters are fit from fit_type (may be centered
;   NFW, single offset NFW, w/ or w/o PS term)
; if center, refcen, and groupFile are given, plot model NFW offset by the distribution of distances between centers
; if fit_type2 and p_mean2 are given, plot two models for comparison

defsysv,'!Omega_m', exists=exists
if NOT exists THEN define_cosmo

;-------------------------------------------------------------------------
; Data Struct
;-------------------------------------------------------------------------
full_str=mrdfits(lensing_infile,1)
sel_str=where(full_str.e1_num GE 10)

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
   x    = full_str.plot_radius_kpc[sel_str]*rnfw ; this will break since rnfw no longer defined here, but stackx is currently not used
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

   get_ds_model,fit_type,p_mean,full_str,x_mpc,ps_term=ps_term,nfw_term=nfw_term,$
                center=center,refcen=refcen,groupFile=groupFile,nfw_off=nfw_off,use_m200=use_m200,$
                mnfw=mnfw,conc=conc,rnfw=rnfw
   ; Baryonic point source term
   if(fit_type[0] NE 0) then oplot,x_mpc,ps_term,color=!red,linestyle=1

   ; NFW term
   oplot,x_mpc,nfw_term,color=!green,linestyle=2

   ; NFW term with offset distribution
   if(keyword_set(center) AND keyword_set(refcen)) then begin
      oplot,x_mpc,nfw_off,color=!orange,linestyle=3
   endif

   ; Sum of terms (currently just NFW + point source if point source is included in model)
   if(fit_type[0] NE 0) then begin
      if(keyword_set(center) AND keyword_set(refcen)) then tot = ps_term + nfw_off $
      else tot = ps_term + nfw_term
      oplot,x_mpc,tot,color=!blue
   endif

   chisq=get_ds_chisq(fit_type,p_mean,full_str,x,y,yerr,center=center,refcen=refcen,groupFile=groupFile,dof=dof,use_m200=use_m200)

   if(keyword_set(fit_type2) AND keyword_set(p_mean2)) then begin
      ;-------------------------------------------------------------------------
      ; PLOT 2ND MODEL
      ;-------------------------------------------------------------------------
      get_ds_model,fit_type2,p_mean2,full_str,x_mpc,ps_term=ps_term,nfw_term=nfw_term,$
                   center=center,refcen=refcen,groupFile=groupFile,nfw_off=nfw_off,use_m200=use_m200,$
                   mnfw=mnfw2,conc=conc2,rnfw=rnfw2

      ; Baryonic point source term
      if(fit_type2[0] NE 0) then oplot,x_mpc,ps_term,color=!red,linestyle=1

      ; NFW term
      oplot,x_mpc,nfw_term,color=!green,linestyle=2

      ; NFW term with offset distribution
      if(keyword_set(center) AND keyword_set(refcen)) then begin
         oplot,x_mpc,nfw_off,color=!orange,linestyle=3
      endif

      ; Sum of terms (currently just NFW + point source if point source is included in model)
      if(fit_type2[0] NE 0) then begin
         if(keyword_set(center) AND keyword_set(refcen)) then tot = ps_term + nfw_off $
         else tot = ps_term + nfw_term
         oplot,x_mpc,tot,color=!blue
      endif

      chisq2=get_ds_chisq(fit_type2,p_mean2,full_str,x,y,yerr,center=center,refcen=refcen,groupFile=groupFile,dof=dof2,use_m200=use_m200)
   endif

   ;-------------------------------------------------------------------------
   ; Replot points
   ;-------------------------------------------------------------------------
   oploterror,x,y,yerr,psym=8,symsize=1.0,color=!black

endif

;-------------------------------------------------------------------------
; Legend
;-------------------------------------------------------------------------
nlens    = textoidl('N_{Lens}:')+string(full_str.lens,format="(I)")
z        = 'Redshift:'+string(full_str.z_lens,format="(f10.2)")
if(keyword_set(models)) then begin
   if(NOT(keyword_set(fit_type2) AND keyword_set(p_mean2))) then begin
      chisq_str = textoidl('\chi^2:')+string(chisq,format='(f10.2)')
      dof_str = 'd.o.f.:'+string(dof,format='(I)')
      if keyword_set(use_m200) then begin
         m        = textoidl('log_{10}(M_{200}):')+string(mnfw,format="(f10.2)")
         r        = textoidl('R_{200}:')+string(rnfw,format="(f10.2)")
      endif else begin
         m        = textoidl('log(M_{vir}):')+string(mnfw,format="(f10.2)")
         r        = textoidl('R_{vir}:')+string(rnfw,format="(f10.2)")
      endelse
      c        = 'Concentration:'+string(Conc,format="(f10.2)")
      items=[nlens,z,m,chisq_str,dof_str]
   endif else begin
      chisq_str = textoidl('\chi^2:')+string(chisq,format='(f6.2)')+','+string(chisq2,format='(f6.2)')
      dof_str = 'd.o.f.:'+string(dof,format='(I)')+','+string(dof2,format='(I)')
      if keyword_set(use_m200) then begin
         m        = textoidl('log_{10}(M_{200}):')+string(mnfw,format="(f6.2)")+','+string(mnfw2,format='(f6.2)')
         r        = textoidl('R_{200}:')+string(rnfw,format="(f6.2)")+','+string(rnfw2,format='(f6.2)')
      endif else begin
         m        = textoidl('log(M_{vir}):')+string(mnfw,format="(f6.2)")+','+string(mnfw2,format='(f6.2)')
         r        = textoidl('R_{vir}:')+string(rnfw,format="(f6.2)")+','+string(rnfw2,format='(f6.2)')
      endelse
      c        = 'Concentration:'+string(Conc,format="(f6.2)")+','+string(conc2,format='(f6.2)')
      roff_str='Offset: 0,'+string(p_mean2[1],format='(F6.2)')
      items=[nlens,z,m,roff_str,chisq_str,dof_str]
   endelse
endif else items=[nlens,z]
legend,items,/right,linestyle=-99,box=0,spacing=1.5
    
ps_close

return
end
