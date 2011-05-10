pro make_plot_offset, zarray, str, full_str, p_mean, p_sigma,title,$
                      xr=xr,$
                      yr=yr,$
                      tit=tit,$
                      xnoedge=xnoedge,$
                      ynoedge=ynoedge,$
                      mcmc=mcmc,$
                      groups=groups,$
                      center=center,$
                      refcen=refcen,$
                      xlog=xlog,$
                      ylog=ylog

defsysv,'!Omega_m', exists=exists
if NOT exists THEN define_cosmo

;---------------------------------------------
;    Common Block
;---------------------------------------------

common fit_options, q_c, lens_redshift, fit_type, lens_m_sun, log_sm, use_m200, neg_points, pos_points,str2,str3, use_group, use_maccio, xmar, ymar, xchars, ychars, no_title, ws_corr, lz_mean, sx

;-------------------------------------------------------------------------
; PLOTTING STUFF
;-------------------------------------------------------------------------

; FILLED CIRCLE  
A = FINDGEN(17) * (!PI*2/16.)   
USERSYM, COS(A), SIN(A), /FILL
simpctable

if n_elements(xtit) eq 0 then begin
    xtit='Physical transverse distance,   R [ Mpc  h!D72!N!E-1!N ]'
endif

if n_elements(ytit) eq 0 then begin
    sun=sunsymbol()
    dunit='  [  h!D72!N   M'+sun+' !N  pc!E -2!N  ]'   
    ytit='!MD!X!MS!X  '+dunit
endif

if (no_title eq 1) then begin
    xtit=''
    ytit=''
endif

if(keyword_set(ylog)) then ytickf='loglabels' else ytickf=''
if(keyword_set(xlog)) then xtickf='loglabels' else xtickf=''

xst=1
yst=1

if keyword_set(xnoedge) then begin
    xtit=''
    xtickf='nolabels'
endif

if keyword_set(ynoedge) then begin
    ytit=''
    ytickf='nolabels'
endif

;-------------------------------------------------------------------------
; Str
;-------------------------------------------------------------------------

if (sx gt 0) then begin
    ; Calculate r200
    if (use_m200 eq 1) then begin
        overdensity = get_overdensity(lens_redshift, /r200)
    endif else begin
        overdensity = density_contrast(lens_redshift) ; virial
    endelse

    rho_crit = critical_density(lens_redshift)
    factor   = alog10(double((4.0/3.0)*!Pi*overdensity*rho_crit))
    r_log    = (1.0/3.0)*(p_mean[0]-factor)
    rnfw     = 10.0^(r_log)

    x    = str3.meanr*rnfw
    y    = str3.sigma
    yerr = str3.sigmaerr

    if (neg_points eq 1) then begin
        str2.meanr = str2.meanr*rnfw*1e3
    endif

    ; now set back to 0
    sx = 0
endif else begin
    x    = str3.meanr/1e3
    y    = str3.sigma
    yerr = str3.sigmaerr
endelse

;-------------------------------------------------------------------------
; Variables
;-------------------------------------------------------------------------

; PP = [ M0, Rnfw, Conc, alpha, bias, mnfw]
PP = ds_model(x,p_mean,/return_all)

M0      = PP[0]
Rnfw    = PP[1] 
Conc    = PP[2] 
alpha   = PP[3] 
bias    = PP[4]
mnfw    = PP[5]  ; in log units
m_sigma = PP[6]

zl      = lens_redshift

;-------------------------------------------------------------------------
; PLOT POINTS
;-------------------------------------------------------------------------

if (xmar[0] eq -99 and ymar[0] eq -99) then begin
    ploterror,x,y,yerr,xlog=xlog,ylog=ylog,title=title,yr=yr,xr=xr,xtickf=xtickf,ytickf=ytickf,$
              tit=tit,xtit=xtit,ytit=ytit,/xst,/yst,psym=3,xchars=xchars,ychars=ychars
endif else begin
    if (xmar[0] ne -99) then begin
        if (ymar[0] ne -99) then begin 
            ploterror,x,y,yerr,xlog=xlog,ylog=ylog,title=title,yr=yr,xr=xr,xtickf=xtickf,ytickf=ytickf,$
                      tit=tit,xtit=xtit,ytit=ytit,/xst,/yst,psym=3,xmargin=xmar,ymargin=ymar,xchars=xchars,ychars=ychars 
        endif else begin
            ploterror,x,y,yerr,xlog=xlog,ylog=ylog,title=title,yr=yr,xr=xr,xtickf=xtickf,ytickf=ytickf,$
                      tit=tit,xtit=xtit,ytit=ytit,/xst,/yst,psym=3,xmargin=xmar,xchars=xchars,ychars=ychars 
        endelse
    endif else begin
        ploterror,x,y,yerr,xlog=xlog,ylog=ylog,title=title,yr=yr,xr=xr,xtickf=xtickf,ytickf=ytickf,$
                  tit=tit,xtit=xtit,ytit=ytit,/xst,/yst,psym=3,xmargin=xmar,ymargin=ymar,xchars=xchars,ychars=ychars
    endelse
endelse

;oploterror,(x+0.05),y,str.sigmaerr_jack,color=!orange,errcolor=!orange

;-------------------------------------------------------------------------
; PLOT DIFFERENCE IN C
;-------------------------------------------------------------------------

;z1 = zarray[0]
;P1  = [p_mean,z1]
;ds2 = ds_model(x,P1,/nfw)
;oplot, x, ds2, color=!grey60,linestyle=1

;z1 = zarray[1]
;P1  = [p_mean,z1]
;ds2 = ds_model(x,P1,/nfw)
;oplot, x, ds2, color=!grey60,linestyle=1

;-------------------------------------------------------------------------
; PLOT THE MODEL
;-------------------------------------------------------------------------

;x_mpc = (findgen(600)+1)/200.0    ; 1 to 3
;x_mpc = findgen(60)/60*(2-0.01)+0.01
nxMpc = 50
if(keyword_set(xlog)) then begin
   xbuffer=1.3
   x_mpc = 10.^(findgen(nxMpc)/(nxMpc-1)*alog10((xr[1]*xbuffer)/(xr[0]/xbuffer)))*x[0]/xbuffer
;   sel = where(x_mpc ge 0.01)
;   x_mpc = x_mpc(sel) 
endif else x_mpc = findgen(nxMpc)/(nxMpc-1) * (xr[1]-xr[0]) + xr[0]

;ds2 = ds_model(x_mpc,p_mean)
;oplot, x_mpc, ds2, color=!blue

nfw = ds_model(x_mpc,p_mean,/nfw)
oplot,x_mpc,nfw,color=!blue,linestyle=2
nfw_off = ds_model(x_mpc,p_mean,center=center,refcen=refcen,/offset,/nfw)
oplot,x_mpc,nfw_off,color=!red,linestyle=0


; NEED TO CHECK THIS HERE .....
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
; Plot negative points
;-------------------------------------------------------------------------

if (neg_points eq 1) then begin
    x2    = str2.meanr/1e3
    if(keyword_set(ylog)) then begin
       y2    = abs(str2.sigma)
       yerr2 = abs(str2.sigmaerr)
       psym=6
       color=!grey30
       errstyle=1
       symsize=0.8
    endif else begin
       y2    = str2.sigma
       yerr2 = str2.sigmaerr
       psym=8
       symsize=1.0
       color=!black
       errstyle=0
    endelse
;    oploterror,x2,y2,yerr2,psym=6,symsize=0.8,color=!grey25,errstyle=1
;    oploterror,x2,y2,yerr2,psym=6,symsize=0.8,color=!grey30,errstyle=1
    oploterror,x2,y2,yerr2,psym=psym,symsize=symsize,color=color,errstyle=errstyle
endif

;-------------------------------------------------------------------------
; Arrow
;-------------------------------------------------------------------------

;x_hl = (full_str.half_light)*10.0/1e3
;arrow,x_hl,1.5,x_hl,0.5 ,/data,/solid,thick=3,color=!orange
;xyouts,(x_hl-0.015),3.2,'half light',color=!orange
;xyouts,(x_hl-0.015),2.0,'radius x 10',color=!orange
;arrow,x_hl,3,x_hl,1.2 ,/data,/solid,thick=3,color=!orange
;xyouts,(x_hl-0.005),3.5,'10 x Rhl',color=!orange

;-------------------------------------------------------------------------
; Star formation efficiency
; Mandelbaum sect 3.1
;-------------------------------------------------------------------------

stellar    = alog10(lens_m_sun*1e12)   ; h^-1 Msun
virialmass = mnfw
 
eta = (stellar/mnfw)*(!Omega_m/!Omega_b)

nlens    = 'N!DLENS!N:'+string(full_str.lens,format="(I)")
if keyword_set(use_m200) then begin
    m        = 'log!D10!N(M!D200!N):'+string(virialmass,format="(f10.2)")
    r        = 'R200:'+string(rnfw,format="(f10.2)")
endif else begin
    m        = 'log(M!Dvir!N):'+string(virialmass,format="(f10.2)")
    r        = 'Rvir:'+string(rnfw,format="(f10.2)")
endelse
c        = 'Concentration:'+string(Conc,format="(f10.2)")
z        = 'Redshift:'+string(lens_redshift,format="(f10.2)")
msun     = 'log!D10!N(M!Dstellar!N):'+string(stellar,format="(f10.2)") 
eta      = 'SF eff:'+string(eta,format="(f10.1)") 
bias     = 'Bias:'+string(bias,format="(f10.2)") 
alpha    = 'Alpha:'+string(alpha,format="(f10.2)") 
m_sigma  = 'm_sigma:'+string(m_sigma,format="(f10.1)") 

;box = flag to include/omit box around the legend (D=include)
;		  outline_color = color of box outline (D = !P.color)


if keyword_set(groups) then begin
  items=[nlens,z,m,c]
  legend,items,/right,linestyle=[-99,-99,-99,-99],box=0,spacing=1.5
endif else begin
  ;items=[nlens,z,msun,m,c,bias,alpha]
  ;legend,items,/right,linestyle=[-99,-99,-99,-99,-99,-99,-99]
  items=[nlens,z,msun,m,c]
  legend,items,/right,linestyle=[-99,-99,-99,-99,-99],box=0,spacing=1.5
endelse
    
return
end

