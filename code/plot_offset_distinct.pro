pro plot_offset_distinct, file, ps_file, mnfw,$
                          center=center,$
                          refcen=refcen,$
                          stackx=stackx,$
                          use_m200=use_m200,$
                          use_maccio=use_maccio

defsysv,'!Omega_m', exists=exists
if NOT exists THEN define_cosmo

full_str=mrdfits(file,1)

;-------------------------------------------------------------------------
; PLOTTING STUFF
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
if(keyword_set(ylog)) then yr = [0.5,3000] else yr = [-200,500]
if(keyword_set(ylog)) then ytickf='loglabels' else ytickf=''
if(keyword_set(xlog)) then xtickf='loglabels' else xtickf=''

;-------------------------------------------------------------------------
; Str
;-------------------------------------------------------------------------
sel_str=where(full_str.e1_num GE 10)

if (keyword_set(stackx)) then begin
    ; Calculate r200
    if (keyword_set(use_m200)) then begin
        overdensity = get_overdensity(full_str.z_lens, /r200)
    endif else begin
        overdensity = density_contrast(full_str.z_lens) ; virial
    endelse

    rho_crit = critical_density(full_str.z_lens)
    factor   = alog10(double((4.0/3.0)*!Pi*overdensity*rho_crit))
    r_log    = (1.0/3.0)*(mnfw-factor)
    rnfw     = 10.0^(r_log)

    x    = full_str.plot_radius_kpc[sel_str]*rnfw
    y    = full_str.we1_mean[sel_str]
    yerr = full_str.we1_error[sel_str]

endif else begin
    x    = full_str.plot_radius_kpc[sel_str]/1e3
    y    = full_str.we1_mean[sel_str]
    yerr = full_str.we1_error[sel_str]
endelse

;-------------------------------------------------------------------------
; Variables
;-------------------------------------------------------------------------

; Assumes 1-parameter fit - this avoids using ds_model and common block
zl = full_str.z_lens
if (use_m200 eq 1) then begin
   if (keyword_set(use_maccio)) then begin
      conc = get_conc(Mnfw,zl,/use200,/maccio) ; Maccio
   endif else begin
      conc = get_conc(Mnfw,zl,/use200) ; Zhao
   endelse
endif else begin
   conc = get_conc(Mnfw,zl,/virial)
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

; PP = [ M0, Rnfw, Conc, alpha, bias, mnfw]
;PP = ds_model(x,p_mean,/return_all,zl=zl)

;M0      = PP[0]
;Rnfw    = PP[1] 
;Conc    = PP[2] 
;alpha   = PP[3] 
;bias    = PP[4]
;mnfw    = PP[5]  ; in log units
;m_sigma = PP[6]


;-------------------------------------------------------------------------
; PLOT POINTS
;-------------------------------------------------------------------------
ploterror,x,y,yerr,xlog=xlog,ylog=ylog,yr=yr,xr=xr,xtickf=xtickf,ytickf=ytickf,$
          xtitle=xtitle,ytitle=ytitle,/xst,/yst,psym=3

;-------------------------------------------------------------------------
; PLOT THE MODEL
;-------------------------------------------------------------------------
nxMpc = 50
if(keyword_set(xlog)) then begin
   xbuffer=1.3
   x_mpc = 10.^(findgen(nxMpc)/(nxMpc-1)*alog10((xr[1]*xbuffer)/(xr[0]/xbuffer)))*xr[0]/xbuffer
endif else x_mpc = findgen(nxMpc)/(nxMpc-1) * (xr[1]-xr[0]) + xr[0]

nfw = nfw_ds(x_mpc,[rnfw,conc],zl,r200=keyword_set(use_m200))
nfw_off = nfw_ds_offset(x_mpc,[rnfw,conc],zl,r200=keyword_set(use_m200),center=center,refcen=refcen)
oplot,x_mpc,nfw_off,color=!red,linestyle=0
oplot,x_mpc,nfw,color=!blue,linestyle=2

;nfw = ds_model(x_mpc,p_mean,/nfw,zl=zl)
;oplot,x_mpc,nfw,color=!blue,linestyle=2
;nfw_off = ds_model(x_mpc,p_mean,center=center,refcen=refcen,/offset,/nfw,zl=zl)
;oplot,x_mpc,nfw_off,color=!red,linestyle=0

;nfw = ds_model(x_mpc,p_mean,/ds_z0)
;oplot,x_mpc,nfw,color=!grey20,linestyle=2

;point_s = ds_model(x_mpc,p_mean,/point_s)
;oplot,x_mpc,point_s,color=!red,linestyle=3

;two_halo = ds_model(x_mpc,p_mean,/two_halo)
;oplot,x_mpc,two_halo,color=!grey30,linestyle=4

;weak_shear = ds_model(x_mpc,p_mean,/weak_nfw)
;oplot,x_mpc,weak_shear,color=!orange,linestyle=4

;group = ds_model(x_mpc,p_mean,/group)
;oplot,x_mpc,group,color=!magenta

; TEMPORARY: ADDED THIS HERE FOR MATTS CENTER STUDY
;p_mean_mmgg_scale = [13.446676]
;nfw_mmgg_scale = ds_model(x_mpc,p_mean_mmgg_scale,/nfw)
;oplot,x_mpc,nfw_mmgg_scale,color=!blue,linestyle=0

; IF WE FIT M0.....
; Plot the real M0
;pointmass    = (10^(log_sm))/1e12   ; 10^12 h^1 Msun
;p_mean2      = p_mean
;p_mean2[4]   = pointmass
;PP2          = [p_mean2,zl]
;point_s      = ds_model(x_mpc,PP2,/point_s)
;oplot,x_mpc,point_s,color=!red, linestyle=2

;-------------------------------------------------------------------------
; Replot points
;-------------------------------------------------------------------------
oploterror,x,y,yerr,psym=8,symsize=1.0,color=!black

;-------------------------------------------------------------------------
; Legend
;-------------------------------------------------------------------------
nlens    = textoidl('N_{Lens}:')+string(full_str.lens,format="(I)")
if keyword_set(use_m200) then begin
    m        = textoidl('log_{10}(M_{200}):')+string(mnfw,format="(f10.2)")
    r        = textoidl('R_{200}:')+string(rnfw,format="(f10.2)")
endif else begin
    m        = textoidl('log(M_{vir}):')+string(mnfw,format="(f10.2)")
    r        = textoidl('R_{vir}:')+string(rnfw,format="(f10.2)")
endelse
c        = 'Concentration:'+string(Conc,format="(f10.2)")
z        = 'Redshift:'+string(full_str.z_lens,format="(f10.2)")
;msun     = textoidl('log_{10}(M_{stellar}):')+string(stellar,format="(f10.2)") 
;eta      = 'SF eff:'+string(eta,format="(f10.1)") 
;bias     = 'Bias:'+string(bias,format="(f10.2)") 
;alpha    = 'Alpha:'+string(alpha,format="(f10.2)") 
;m_sigma  = 'm_sigma:'+string(m_sigma,format="(f10.1)") 

;box = flag to include/omit box around the legend (D=include)
;		  outline_color = color of box outline (D = !P.color)

items=[nlens,z,m,c]
legend,items,/right,linestyle=[-99,-99,-99,-99],box=0,spacing=1.5
    
ps_close

return
end
