function find_contour_levels,res
; get contour levels that enclose 68%, 95%, 99% of mock data points

sel68=where(total(res[reverse(sort(res))],/cumulative) GT 0.68*total(res),n68)
sel95=where(total(res[reverse(sort(res))],/cumulative) GT 0.95*total(res),n95)
sel99=where(total(res[reverse(sort(res))],/cumulative) GT 0.99*total(res),n99)

if(n68 GT 0) then level1=res[(reverse(sort(res)))[min(sel68)]]
if(n95 GT 0) then level2=res[(reverse(sort(res)))[min(sel95)]]
if(n99 GT 0) then level3=res[(reverse(sort(res)))[min(sel99)]]

if(n_elements(level1) GT 0 $
   AND n_elements(level2) GT 0 $
   AND n_elements(level3) GT 0) then begin
      if(level1 GT level2 AND level2 GT level3) then $
         return,[level3,level2,level1] $
      else return,-1
endif else return,-1

end

pro oplot_contours, x,y,xBin,yBin,xRange,yRange,color,g_smooth=g_smooth
; overplot contours for mock distribution onto 2-d 
res=hist_2d(x,y,min1=min(xRange),min2=min(yRange),max1=max(xRange),max2=max(yRange),bin1=xBin,bin2=yBin)

if n_elements(g_smooth) gt 0 then begin	;make a gaussian and convolve with it
   make_gaussian,gauss,size=[2.0*g_smooth,2.0*g_smooth],$
                 fwhm=g_smooth,counts=1.0
   res=convol(res,gauss,/edge_tr)
endif

dims=size(res,/dim) ;get dimensions of the 2d histogram
xs=findgen(dims[0])*xBin+min(xRange)
ys=findgen(dims[1])*yBin+min(yRange)
levels=find_contour_levels(res)
if(n_elements(levels) GT 1) then $ 
   contour,res,xs,ys,color=color,/over,levels=levels[1:2],/fill,c_color=[!medblue,!darkblue]
; else plot no contours
end

pro ds_cov_plots, chainFile, fit_type, plotFile, hist=hist, burnin=burnin

; copied from /Users/alexie/idl/MCMC/dsmodel/ds_cov_plots.pro with changes
;-----------------------------------------------------------------------------
; Adapted from Dave : /Users/alexie/idl/dave_idl/inversion/mcmc_cov_plots3.pro
;-----------------------------------------------------------------------------

mcmc=obj_new('mcmc')
pars=mcmc->read_trials(chainFile)
if(n_elements(burnin) EQ 0) then burnin=500
pars=pars[*,burnin:*]

sz=size(pars,/dim)   ; pars = chain
npars=sz[0]              ; n=number of pars

; -- Set up plot Labels -------------------------------------------------
; 0  baryonic mass
; 1  R_vir : NFW virial radius
; 2  C     : NFW concentration
; 3  alpha : fraction
; 4  Bias
; 5  m_sigma : dispersion in central mass
;-------------------------------------------------------------------------

sun=sunsymbol()
titles = strarr(npars)
ranges=fltarr(2,npars)
tickv=fltarr(3,npars)
nticks=n_elements(tickv[*,0])-1
minor=intarr(npars)

k = 0
; Baryonic term
if (fit_type[0] eq 1) then begin
    titles[k] = 'log(M_{cen}/M'+sun+')'
    ranges[*,k]=[9.8,11.8]
    tickv[*,k]=[10.0,10.8,11.6]
    minor[k]=8
    k=k+1  
endif
; Mnfw
if (fit_type[1] eq 1) then begin 
    titles[k] = 'log(M_{200c}/M'+sun+')'
    ranges[*,k]=[13.,13.8]
    tickv[*,k]=[13.1,13.4,13.7]
    minor[k]=3
    k=k+1  
endif
; Conc
if ((fit_type[2] eq 1)) then begin
    titles[k] = 'c_{200c}'
    ranges[*,k]=[1,12]
    tickv[*,k]=[2.,6.,10.]
    minor[k]=4
    k=k+1  
endif
; Alpha
if (fit_type[3] eq 1) then begin 
    titles[k] = 'Alpha'
    k=k+1  
endif
; Bias
if (fit_type[4] eq 1) then begin 
    titles[k] = 'Bias'
    k=k+1  
endif
; Dispersion in central mass
if (fit_type[5] eq 1) then begin 
    titles[k] = 'Dispersion'
    k=k+1  
endif
; Offset radius
if (fit_type[6] eq 1) then begin 
    titles[k] = 'R_{off} (kpc)'
    ranges[*,k]=[0,100]
    tickv[*,k]=[0.,50.,100.]
    minor[k]=5
    pars[k,*]*=1000. ; convert from Mpc to kpc
endif

titles=textoidl(titles)

;-------------------------------------------------------------------------
; Marginal Posterior distribution for parameter pairs
;------------------------------------------------------------------------- 
    
!p.multi=[0,npars,npars]
!p.thick=3
!x.thick=3
!y.thick=3
!p.charsize=2
!p.charthick=1.2
!p.font=0

; position setup
lmarg=0.1
rmarg=0.05
bmarg=0.1
tmarg=0.05
xspace=0.0
yspace=0.0
xwidth=(1.-lmarg-rmarg-xspace*(npars-2))/(npars-1)
ywidth=(1.-bmarg-tmarg-yspace*(npars-2))/(npars-1)

set_plot,'ps'
simpctable
device,filename=plotFile,/encapsul,/color,/helvetica,xsize=8,ysize=8,/inches

for ii=0, npars-1 do begin            ; Go through params
   for jj=ii+1, npars-1 do begin      ; Go through params
      

      blank=replicate(' ',nticks+1)
      fill=replicate('',nticks+1)
      if(jj EQ npars-1) then begin
         xtit=titles[ii]
         xtickn=fill
      endif else begin
         xtit=''
         xtickn=blank
      endelse
      if(ii EQ 0) then begin
         ytit=titles[jj]
         ytickn=fill
      endif else begin
         ytit=''
         ytickn=blank
      endelse
      
      ; position
      x1=lmarg+ii*(xwidth+xspace)
      x2=x1+xwidth
      y1=bmarg+(npars-jj-1)*(ywidth+yspace)
      y2=y1+ywidth

      xw=sigma(pars[ii,*])
      yw=sigma(pars[jj,*])
      nn=n_elements(pars[ii,*])
      xp=pars[ii,*]
      yp=pars[jj,*]
      xran=randomn(seed,nn)*xw/5.0
      yran=randomn(seed,nn)*yw/5.0
                
      xr=ranges[*,ii]
      yr=ranges[*,jj]
      nBins=50.
      xbin=(xr[1]-xr[0])/nBins
      ybin=(yr[1]-yr[0])/nBins

      plot,/nodata,xr,yr,xst=1+4,yst=1+4,position=[x1,y1,x2,y2]
      oplot_contours,xp,yp,xbin,ybin,xr,yr,!black ;,g_smooth=sm
      axis,xaxis=0,xst=1,xtickn=xtickn,xtickv=tickv[*,ii],xticks=nticks,xminor=minor[ii],xtit=xtit
      axis,xaxis=1,xst=1,xtickn=blank,xtickv=tickv[*,ii],xticks=nticks,xminor=minor[ii]
      axis,yaxis=0,yst=1,ytickn=ytickn,ytickv=tickv[*,jj],yticks=nticks,yminor=minor[jj],ytit=ytit
      axis,yaxis=1,yst=1,ytickn=blank,ytickv=tickv[*,jj],yticks=nticks,yminor=minor[jj]
   endfor
endfor
    
device,/close

return
end



