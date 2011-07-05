pro test_offset_range

; try plotting lensing signal from MMGG_scale along with some offset
; models to see what range of offsets are still feasible given the
; data.

; some code copied from plot_lensing_results.pro

str=mrdfits("~/data/cosmos/groups_lensing/outfiles/bin_50_2000_8_emp/center_mmgg_scale.fits",1)

zl=str.z_lens
mnfw=13.37
overdensity = get_overdensity( zl, /r200)
rho_crit = critical_density(zl)
factor   = alog10(double((4.0/3.0)*!Pi*overdensity*rho_crit))
r_log    = (1.0/3.0)*(Mnfw-factor)
rnfw     = 10.0^(r_log)
conc = get_conc(Mnfw,zl,/use200) ; Zhao

x=str.plot_radius_kpc*rnfw
y=str.we1_mean
yerr=str.we1_error

!p.font=0
!p.thick=3
!x.thick=3
!y.thick=3
!p.charsize=1.2
!p.charthick=1.2
set_plot,'ps'
device,filename='test_offset_range.eps',/encapsul,/color

xst=1
yst=1
xlog=1
ylog=0
if(keyword_set(xlog)) then xr = [0.02,2] else xr=[0,1.5]
if(keyword_set(ylog)) then yr = [0.5,3000] else yr = [-200,400]
if(keyword_set(ylog)) then ytickf='loglabels' else ytickf=''
if(keyword_set(xlog)) then xtickf='loglabels' else xtickf=''

plotcircle
simpctable

xtitle=textoidl('Physical transverse distance,  R (h_{72}^{-1} Mpc)')
bar='!S!A=!R!N'
ytitle=textoidl('\Delta\Sigma = ')+bar+textoidl('\Sigma(<R) - \Sigma(R)  (h_{72} M')+sunsymbol()+textoidl(' pc^{-2})')

ploterror,x,y,yerr,xlog=xlog,ylog=ylog,yr=yr,xr=xr,xtickf=xtickf,ytickf=ytickf,$
          xtitle=xtitle,ytitle=ytitle,/xst,/yst,psym=8,symsize=1.0,color=!black


; overplot models
   nxMpc = 50
   if(keyword_set(xlog)) then begin
      xbuffer=1.3
      x_mpc = 10.^(findgen(nxMpc)/(nxMpc-1)*alog10((xr[1]*xbuffer)/(xr[0]/xbuffer)))*xr[0]/xbuffer
   endif else x_mpc = findgen(nxMpc)/(nxMpc-1) * (xr[1]-xr[0]) + xr[0]
   
   ; NFW term
   nfw = nfw_ds(x_mpc,[rnfw,conc],zl,/r200)
   oplot,x_mpc,nfw,color=!green,linestyle=2

   roff=[0.01,0.02,0.03,0.04,0.05] ; Mpc
   colors=[!red,!darkorange,!orange,!darkyellow,!yellow]
   for rr=0,n_elements(roff)-1 do begin

      sigmaR = nfw_sigma_offset(x_mpc,roff[rr],[rnfw,conc],zl,/r200)
      sigmaMean=dblarr(n_elements(x_mpc))
      
   ; \bar{\Sigma}(<R|R_{off})
      for i=0,n_elements(x_mpc)-1 do begin
         npts=300               ; slow, but needed for decent convergence at small r
         dr=double(x_mpc[i])/npts
         rInterior = dindgen(npts)*dr
         sigmaInterior = nfw_sigma_offset(rInterior,roff[rr],[rnfw,conc],zl,/r200)
         integrand=rInterior*sigmaInterior ; sigmaInterior goes to infinity at r=0, so remove
         sigmaMean[i] = 2./x_mpc[i]^2 * total(integrand[where(finite(integrand))]*dr)
      endfor

      deltaSigma=sigmaMean-sigmaR

      oplot,x_mpc,deltaSigma,color=colors[rr],linestyle=3
   endfor
device,/close

stop
end
