pro plot_full_stacks,cenText,lensFileArr,plotFile,massMean,fitTypeAll,stackx=stackx,use_m200=use_m200,test=test

; plot the full lensing stacks with model fits for different centers
; in separate panels
; builds on plot_lensing_results.pro but combines all plots for different
; individual centers into one big multi-panel plot
  
; assumes mass is the only fit parameter and adds central galaxy point
; source according to which center is used

; MRG 5 July 2011

if(keyword_set(test)) then begin
   ; Set paths for output files
   dirName='bin_50_1500_7_emp'
   fileDir='~/data/cosmos/groups_lensing/outfiles/'+dirName+'/'
   plotDir='~/data/cosmos/groups_lensing/plots/'+dirName+'/'

   cenNames=['mmgg_scale','mmgg_r200','mlgg_scale','mlgg_r200','cl','cm','xray','cn']
   cenText=textoidl(['MMGG_{scale}','MMGG_{R200}','MLGG_{scale}','MLGG_{R200}','CL','CM','X-ray','CN'])
   ptSrc=[2,2,2,2,0,0,0,0]      ; for fit_t

   lensFileArr=strcompress(fileDir+'center_'+cenNames+'.fits',/remove_all)
   plotFile=plotDir+'full_stacks.eps'

   ; create 2d array to save fit types, 1 row for each center
   fit_t = [$
           2,$                  ; 0  M0    : baryonic mass
           1,$                  ; 1  R_vir : NFW virial mass
           2,$                  ; 2  C     : NFW concentration
           0,$                  ; 5  alpha : fraction
           0,$                  ; 6  bias
           0 ]                  ; 6  m_sigma
   fitTypeAll=rebin(fit_t,n_elements(fit_t),n_elements(cenNames))
   fitTypeAll[0,*]=ptSrc

   ; arrays to record mean and stddev mass from mcmc
   massMean=replicate(13.0,n_elements(cenNames))
endif


nCols=2
nRows=4
!p.multi=[0,nCols,nRows]
!p.font=0
!p.charsize=1.2
!p.charthick=1.2
!p.thick=3
!x.thick=3
!y.thick=3
set_plot,'ps'
simpctable
plotcircle,0.7
sun=sunsymbol()
star=starsymbol()
blank=replicate(' ',60)

;margins and plot dimensions
xstart=0.12
ystart=0.07
dx=(1.-1.2*xstart)/nCols
dy=(1.-1.4*ystart)/nRows

; axes
xst=1
yst=1
xlog=1
ylog=0
if(keyword_set(xlog)) then xr = [0.02,2] else xr=[0,1.5]
if(keyword_set(ylog)) then yr = [0.5,3000] else yr = [-99,200]
;if(keyword_set(ylog)) then ytickf='loglabels' else ytickf=''
;if(keyword_set(xlog)) then xtickf='loglabels' else xtickf=''
xtitle=textoidl('Physical transverse distance,  R (h_{72}^{-1} Mpc)')
bar='!S!A=!R!N'
ytitle=textoidl('\Delta\Sigma = ')+bar+textoidl('\Sigma(<R) - \Sigma(R)  (h_{72} M')+sunsymbol()+textoidl(' pc^{-2})')

device,filename=plotFile,/encapsul,/helvetica,/color,xsize=7,ysize=9,/inches

nCen=n_elements(cenText)
for ii=0,nCen-1 do begin
   ; positions for this plot
   x1=xstart+(ii MOD nCols)*dx
   x2=x1+dx
   y1=ystart+(ii-(ii MOD nCols))/2*dy
   y2=y1+dy

   ; read data file and get parameters
   str=mrdfits(lensFileArr[ii],1)
   conc=get_conc(massMean[ii],str.z_lens,use200=keyword_set(use_m200)) ; Zhao
   overdensity=get_overdensity(str.z_lens,r200=keyword_set(use_m200))
   rho_crit=critical_density(str.z_lens)
   factor=alog10(double((4.0/3.0)*!Pi*overdensity*rho_crit))
   r_log=(1.0/3.0)*(massMean[ii]-factor)
   rnfw=10.0^(r_log)

   ; restrict to points with enough sources
   sel=where(str.e1_num GE 10)
   if(keyword_set(stackx)) then begin
      x    = str.plot_radius_kpc[sel]*rnfw
      y    = str.we1_mean[sel]
      yerr = str.we1_error[sel]
   endif else begin
      x    = str.plot_radius_kpc[sel]/1e3
      y    = str.we1_mean[sel]
      yerr = str.we1_error[sel]
   endelse

   ; Suppress ticknames where necessary
   if(ii+1 GT nCols) then xtickname=blank else xtickname=''
   if((ii MOD nCols) NE 0) then ytickname=blank else ytickname=''

   ; PLOT BLANK PANEL AND AXIS LABELS
   plot,/nodata,xr,yr,position=[x1,y1,x2,y2],xlog=xlog,ylog=ylog,xst=xst,yst=yst,xtitle='',ytitle='',xtickname=xtickname,ytickname=ytickname,xcharsize=1.5,ycharsize=1.5
;   plot,x,y,xr=xr,yr=yr,position=[x1,y1,x2,y2]
   if(ii EQ 0) then begin
      xyouts,xstart+0.5*(nCols*dx),0.25*ystart,xtitle,alignment=0.5,/normal
      xyouts,0.3*xstart,ystart+0.5*(nRows*dy),ytitle,alignment=0.5,orientation=90,/normal
   endif

   ; PLOT THE MODEL
   nxMpc = 50
   if(keyword_set(xlog)) then begin
      xbuffer=1.3
      x_mpc=10.^(findgen(nxMpc)/(nxMpc-1)*alog10((xr[1]*xbuffer)/(xr[0]/xbuffer)))*xr[0]/xbuffer
   endif else x_mpc = findgen(nxMpc)/(nxMpc-1) * (xr[1]-xr[0]) + xr[0]
   
   ; NFW term
   nfw=nfw_ds(x_mpc,[rnfw,conc],str.z_lens,r200=keyword_set(use_m200))

   ; Baryonic point source term
   if(fitTypeAll[0,ii] NE 0) then begin
      ps_term=10^(str.msun_lens)/1.e12/(!pi*x_mpc^2) ; h^-1 Msun, factor of 1e12 to convert to pc^2
      tot=ps_term + nfw
   endif else tot=nfw

   oplot,x_mpc,tot,color=!blue,thick=3
   oplot,x_mpc,nfw,color=!darkgreen,linestyle=2,thick=8
   if(fitTypeAll[0,ii] NE 0) then oplot,x_mpc,ps_term,color=!red,linestyle=1,thick=4


   ; PLOT DATA POINTS
   oploterror,x,y,yerr,psym=8,color=!black
   
   ; CALCULATE CHI^2
   ; currently just NFW + point source if point source is included in model

   ; calculate expected model values at locations of data points
   model_nfw=nfw_ds(x,[rnfw,conc],str.z_lens,r200=keyword_set(use_m200))

   if(fitTypeAll[0,ii] NE 0) then begin
      model_ps=10^(str.msun_lens)/1.e12/(!pi*x^2) ; h^-1 Msun, factor of 1e12 to convert to pc^2
      model_tot=model_ps + model_nfw
   endif else model_tot=model_nfw

   chisq = total((model_tot-y)^2/yerr^2)
   dof = n_elements(x)-n_elements(where(fitTypeAll[*,ii] EQ 1))

   ; LEGEND
   nlens=textoidl('N_{Lens}:')
   zstr='Redshift:'
   if keyword_set(use_m200) then begin
      mstr=textoidl('log(M_{200c}/M'+sun+'):')
      rstr=textoidl('R_{200c} (Mpc):')
   endif else begin
      mstr=textoidl('log(M_{vir}):')
      rstr=textoidl('R_{vir}:')
   endelse
   cstr='Concentration:'
   chisqstr=textoidl('\chi^2:')
   dof_str='d.o.f.:'
   if(fitTypeAll[0,ii] NE 0) then begin
      smstr=textoidl('log(M')+star+'/M'+sun+'):'
   endif else sm=textoidl('log(M')+star+'/M'+sun+'):'+string(0.0,format="(f6.2)")

   yLine=0.1*(yr[1]-yr[0])
   xRight=0.8*xr[1]
   xStrRight=0.33*xr[1]
   fmt="(F6.2)"
   xyouts,xStrRight,yr[1]-1.*yLine,mstr,alignment=0.99,charsize=0.9
   xyouts,xStrRight,yr[1]-2.*yLine,smstr,alignment=0.88,charsize=0.9
   xyouts,xStrRight,yr[1]-3.*yLine,chisqstr,alignment=1.0,charsize=0.9
   xyouts,xRight,yr[1]-1.*yLine,string(massMean[ii],format=fmt),alignment=1,charsize=0.9
   xyouts,xRight,yr[1]-2.*yLine,string(str.msun_lens,format=fmt),alignment=1,charsize=0.9
   xyouts,xRight,yr[1]-3.*yLine,string(chisq,format=fmt),alignment=1,charsize=0.9

   xyouts,xRight,yr[0]+0.5*yLine,cenText[ii],alignment=1

endfor
device,/close

end
