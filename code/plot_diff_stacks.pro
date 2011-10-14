pro plot_diff_stacks,cenNames,refNames,cenText,refText,lensFileArrCen,lensFileArrRef,plotFile,groupFile,stackx=stackx,use_m200=use_m200,test=test

; plot the lensing stacks on different centers vs. mmgg_scale with model fits
; in separate panels, along with histogram of offsets
; builds on plot_lensing_results.pro but combines all plots for different
; individual centers into one big multi-panel plot
  
; assumes mass is the only fit parameter and adds central galaxy point
; source according to which center is used

; MRG 6 July 2011

if(keyword_set(test)) then begin
   ; Set paths for output files
   dirName='bin_50_1500_7_emp'
   fileDir='~/data/cosmos/groups_lensing/outfiles/'+dirName+'/'
   plotDir='~/data/cosmos/groups_lensing/plots/'+dirName+'/'

   cenNames=['mmgg_r200','mlgg_r200','mlgg_scale','cm','cl','cn','xray']
   cenText=textoidl(['MMGG_{R200}','MLGG_{R200}','MLGG_{scale}','CM','CL','CN','X-ray'])
   refNames=replicate('mmgg_scale',n_elements(cenNames)) ; the "good center" to compare with the ones above
   refText=textoidl(replicate('MMGG_{scale}',n_elements(cenNames)))

   lensFileArrCen=strcompress(fileDir+'center_'+cenNames+'_'+refNames+'.fits',/remove_all)
   lensFileArrRef=strcompress(fileDir+'center_'+refNames+'_'+cenNames+'.fits',/remove_all)
   plotFile=plotDir+'diff_stacks_test.eps'
   groupFile='~alexie/Work/GroupCatalogs/cosmos_xgroups_20110209.fits'
endif

; group file needs to be read in for histogram of offsets
group=mrdfits(groupFile,1)
sel=where(group.flag_include EQ 1)
group=group[sel]

nCols=3
nRows=7
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
xstart=0.09
ystart=0.06
ddx=0.01 ; small gap between lensing plots and histogram
dx=(1.-2.*xstart-ddx)/nCols
dy=(1.-1.4*ystart)/nRows

; axes
xst=1
yst=1
xlog=1
ylog=0
yloghist=1
if(keyword_set(xlog)) then xr=[0.02,2] else xr=[0,1.5]
if(keyword_set(ylog)) then begin
   yr=[0.5,3000] 
endif else begin 
   yr=[-120,400]
   ytickv=[0,200,400]
endelse
if(keyword_set(yloghist)) then yrhist=[0.5,150] else yrhist=[0,110]
;if(keyword_set(ylog)) then ytickf='loglabels' else ytickf=''
;if(keyword_set(xlog)) then xtickf='loglabels' else xtickf=''
xtitle=textoidl('Physical transverse distance,  R (h_{72}^{-1} Mpc)')
bar='!S!A=!R!N'
ytitle=textoidl('\Delta\Sigma = ')+bar+textoidl('\Sigma(<R) - \Sigma(R)  (h_{72} M')+sunsymbol()+textoidl(' pc^{-2})')
ytitleHist=textoidl('Distribution of Offsets (N_{groups})')

; legend
yLine=0.12*(yr[1]-yr[0])
xRight=0.8*xr[1]
xStrRight=0.33*xr[1]
fmt="(F6.2)"
lCharSize=0.7

tickCharSize=1.3

device,filename=plotFile,/encapsul,/helvetica,/color,xsize=7,ysize=9,/inches

nCen=n_elements(cenText)
for ii=0,nCen-1 do begin

   ; LEFT COLUMN - MMGG_SCALE
   ; positions for this plot
   x1=xstart
   x2=x1+dx
   y1=ystart+ii*dy
   y2=y1+dy

   ; read data file and get parameters
   str=mrdfits(lensFileArrRef[ii],1)

   fitTypeRef=str.fit_type
   pMeanRef=str.p_mean
   pSigmaRef=str.p_sigma

   ; restrict to points with enough sources
   sel=where(str.e1_num GE 10 AND str.plot_radius_kpc GT 10)
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
   if(ii EQ 0) then begin
      xtickname=''
   endif else begin 
      xtickname=blank
   endelse
   xtit=''
   ytit=''
   ytickname=''

   ; PLOT BLANK PANEL AND AXIS LABELS
   plot,/nodata,xr,yr,position=[x1,y1,x2,y2],xlog=xlog,ylog=ylog,xst=xst,yst=yst,xtitle=xtit,ytitle=ytit,xtickname=xtickname,ytickname=ytickname,xcharsize=tickCharSize,ycharsize=tickCharSize,ytickv=ytickv,yticks=n_elements(ytickv)-1,yminor=4
   if(ii EQ 0) then begin
      xyouts,xstart+0.5*(nCols*dx+ddx),0.25*ystart,xtitle,alignment=0.5,/normal
      xyouts,0.33*xstart,ystart+0.5*(nRows*dy),ytitle,alignment=0.5,orientation=90,/normal
      xyouts,1.73*xstart+nCols*dx+ddx,ystart+0.5*(nRows*dy),ytitleHist,alignment=0.5,orientation=90,/normal
   endif

   ; PLOT THE MODEL
   nxMpc = 50
   if(keyword_set(xlog)) then begin
      xbuffer=1.3
      x_mpc=10.^(findgen(nxMpc)/(nxMpc-1)*alog10((xr[1]*xbuffer)/(xr[0]/xbuffer)))*xr[0]/xbuffer
   endif else x_mpc = findgen(nxMpc)/(nxMpc-1) * (xr[1]-xr[0]) + xr[0]
   
   get_ds_model, fitTypeRef, pMeanRef, str, x_mpc, ps_term=ps_term, nfw_term=nfw_term,$
                 use_m200=use_m200,mnfw=mnfw,conc=conc,rnfw=rnfw  

   ; Sum of terms
   if(fitTypeRef[0] NE 0) then tot=ps_term + nfw_term $
   else tot=nfw_term
   oplot,x_mpc,tot,color=!blue,thick=3

   ; NFW term
   oplot,x_mpc,nfw_term,color=!darkgreen,linestyle=2,thick=8

   ; Baryonic point source term
   if(fitTypeRef[0] NE 0) then oplot,x_mpc,ps_term,color=!red,linestyle=1,thick=6


   ; PLOT DATA POINTS
   oploterror,x,y,yerr,psym=8,color=!black

   ; CALCULATE CHI^2
   chisq=get_ds_chisq(fitTypeRef,pMeanRef,str,x,y,yerr,dof=dof,use_m200=use_m200)

   ; LEGEND
   nstr=textoidl('N_{Lens}:')
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
   if(fitTypeRef[0] NE 0) then begin
      smstr=textoidl('log(M')+star+'/M'+sun+'):'
   endif else sm=textoidl('log(M')+star+'/M'+sun+'):'+string(0.0,format="(f6.2)")

   xyouts,xStrRight,yr[1]-1.*yLine,mstr,alignment=0.99,charsize=lCharSize
   xyouts,xStrRight,yr[1]-2.*yLine,smstr,alignment=0.88,charsize=lCharSize
   xyouts,xStrRight,yr[1]-3.*yLine,chisqstr,alignment=1.0,charsize=lCharSize
   xyouts,xRight,yr[1]-1.*yLine,string(massMeanRef[ii],format=fmt),alignment=1,charsize=lCharSize
   xyouts,xRight,yr[1]-2.*yLine,string(str.msun_lens,format=fmt),alignment=1,charsize=lCharSize
   xyouts,xRight,yr[1]-3.*yLine,string(chisq,format=fmt),alignment=1,charsize=lCharSize

   xyouts,xRight,yr[0]+0.5*yLine,refText[ii],alignment=1,charsize=1.2*lCharSize


   ; MIDDLE COLUMN - OFFSET CENTERS
   ; positions for this plot
   x1=xstart+dx
   x2=x1+dx
   y1=ystart+ii*dy
   y2=y1+dy

   ; read data file and get parameters
   str=mrdfits(lensFileArrCen[ii],1)
   
   ; restrict to points with enough sources
   sel=where(str.e1_num GE 10 AND str.plot_radius_kpc GT 10)
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
   if(ii EQ 0) then begin
      xtickname=''
   endif else begin 
      xtickname=blank
   endelse
   xtit=''
   ytit=''
   ytickname=blank

   ; PLOT BLANK PANEL AND AXIS LABELS
   plot,/nodata,xr,yr,position=[x1,y1,x2,y2],xlog=xlog,ylog=ylog,xst=xst,yst=yst,xtitle=xtit,ytitle=ytit,xtickname=xtickname,ytickname=ytickname,xcharsize=tickCharSize,ycharsize=tickCharSize,ytickv=ytickv,yticks=n_elements(ytickv)-1,yminor=4

   ; PLOT THE MODEL
   nxMpc = 50
   if(keyword_set(xlog)) then begin
      xbuffer=1.3
      x_mpc=10.^(findgen(nxMpc)/(nxMpc-1)*alog10((xr[1]*xbuffer)/(xr[0]/xbuffer)))*xr[0]/xbuffer
   endif else x_mpc = findgen(nxMpc)/(nxMpc-1) * (xr[1]-xr[0]) + xr[0]

                                ; NFW offset model comes from pMeanRef
                                ; with offset distribution from
                                ; groupFile. PS term comes from
                                ; stellar mass saved in str.

   get_ds_model, fitTypeRef, pMeanRef, str, x_mpc, ps_term=ps_term,$
                 center=cenNames[ii],refcen=refNames[ii],groupFile=groupFile,nfw_off=nfw_off, $
                 use_m200=use_m200,mnfw=mnfw,conc=conc,rnfw=rnfw
   
   ; Sum of terms
   if(str.msun_lens GT 0.) then tot=ps_term + nfw_off $
   else tot=nfw_off
   oplot,x_mpc,tot,color=!blue,thick=3

   ; Offset NFW term
   oplot,x_mpc,nfw_off,color=!orange,linestyle=3,thick=8

   ; Baryonic point source term
   if(str.msun_lens GT 0.) then oplot,x_mpc,ps_term,color=!red,linestyle=1,thick=6


   ; PLOT DATA POINTS
   oploterror,x,y,yerr,psym=8,color=!black

   ; CALCULATE CHI^2
;   chisq=get_ds_chisq(fitTypeAllCen[*,ii],massMeanRef[ii],str,x,y,yerr,dof=dof,$
;                      center=cenNames[ii],refcen=refNames[ii],groupFile=groupFile,use_m200=use_m200)

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
   if(str.msun_lens GT 0.) then begin
      smstr=textoidl('log(M')+star+'/M'+sun+'):'
   endif else sm=textoidl('log(M')+star+'/M'+sun+'):'+string(0.0,format="(f6.2)")

   xyouts,xStrRight,yr[1]-1.*yLine,smstr,alignment=0.88,charsize=lCharSize
;   xyouts,xStrRight,yr[1]-2.*yLine,chisqstr,alignment=1.0,charsize=lCharSize
   xyouts,xRight,yr[1]-1.*yLine,string(str.msun_lens,format=fmt),alignment=1,charsize=lCharSize
;   xyouts,xRight,yr[1]-2.*yLine,string(chisq,format=fmt),alignment=1,charsize=lCharSize

   xyouts,xRight,yr[0]+0.5*yLine,cenText[ii],alignment=1,charsize=1.2*lCharSize


   ; RIGHT COLUMN - HISTOGRAM OF OFFSETS
   x1=xstart+2.*dx+ddx
   x2=x1+dx
   y1=ystart+ii*dy
   y2=y1+dy

   get_center_coords,group,cenNames[ii],ra1,dec1,sel1
   get_center_coords,group,refNames[ii],ra2,dec2,sel2
   offset_mpc=distance(ra1,dec1,ra2,dec2)*3600.*group.lensing_r200_mpc/group.lensing_r200_as ; Mpc
   min_offset_mpc=0.050 ; 50 kpc
   bin=0.1
   
   hist=histogram(alog10(offset_mpc),locations=logx,min=alog10(min_offset_mpc),bin=bin)

   ; Suppress ticknames where necessary
   if(ii EQ 0) then begin
      xtickname=''
   endif else begin 
      xtickname=blank
   endelse

   ; PLOT BLANK PANEL AND AXIS LABELS
   plot,/nodata,xr,yrhist,position=[x1,y1,x2,y2],xlog=xlog,ylog=yloghist,xst=1+4,yst=1+4,xtitle=xtit,ytitle=ytit,xtickname=xtickname,ytickname=ytickname

   x=10.^(logx+bin/2.)
   xmin=10.^(logx[0])
   xmax=10.^(logx[n_elements(logx)-1]+bin)
   xpad=[xmin,xmin,x,xmax,xmax]
   ypad=[yrhist[0],hist[0],hist,hist[n_elements(hist)-1],yrhist[0]]
   oplot,xpad,ypad,ps=10,color=!orange

   nMatch=n_elements(where(offset_mpc LT min_offset_mpc))
   xarr=[xr[0],xr[0],min_offset_mpc,min_offset_mpc]
   yarr=[yrhist[0],nMatch,nMatch,yrhist[0]]
   oplot,xarr,yarr,color=!gray,thick=5
   polyfill,xarr,yarr,/line_fill,orientation=30,color=!gray

   ; LEGEND
   nstr=textoidl('N_{Lens}:')
   zstr=textoidl('<z_{Lens}>:')

   xyouts,xStrRight,0.5*yrhist[1],nstr,alignment=1.0,charsize=lCharSize
   xyouts,xStrRight,0.28*yrhist[1],zstr,alignment=1.0,charsize=lCharSize
   xyouts,xRight,0.5*yrhist[1],string(str.lens,format="(I6)"),alignment=1,charsize=lCharSize
   xyouts,xRight,0.28*yrhist[1],string(str.z_lens,format="(F6.2)"),alignment=1,charsize=lCharSize

   axis,xaxis=0,xr=xr,xlog=xlog,xst=xst,xtickname=xtickname,xcharsize=tickCharSize
   axis,xaxis=1,xr=xr,xlog=xlog,xst=xst,xtickname=blank
   axis,yaxis=0,yr=yrhist,ylog=yloghist,yst=yst,ytickname=blank
   axis,yaxis=1,yr=yrhist,ylog=yloghist,yst=yst,ytickname='',ycharsize=tickCharSize
endfor
device,/close

end
