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
   contour,res,xs,ys,color=color,/over,levels=levels[1:2],/fill,c_color=[!red,!blue]
; else plot no contours
end

pro ds_cov_plots, chainFile, fit_type, p_mean, p_sigma, plotFile, hist=hist, burnin=burnin

; copied from /Users/alexie/idl/MCMC/dsmodel/ds_cov_plots.pro with changes
;-----------------------------------------------------------------------------
; Adapted from Dave : /Users/alexie/idl/dave_idl/inversion/mcmc_cov_plots3.pro
;-----------------------------------------------------------------------------

mcmc=obj_new('mcmc')
pars=mcmc->read_trials(chainFile)
if(n_elements(burnin) EQ 0) then burnin=500
pars=pars[*,burnin:*]

;common mcmc_block, x, y, ivar, psigma, npars, Prior_mean, Prior_sigma
;common fit_options, q_c, lens_redshift, fit_type, lens_m_sun

;-- Get robust measures
;mcmc_stats,pars,p_mean,p_sigma,p_mean_rob,p_sigma_rob

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

titles = strarr(npars)

k = 0

; Baryonic term
if (fit_type[0] eq 1) then begin
    titles[k] = 'Mcen'
    k=k+1  
endif
; Rvir
if (fit_type[1] eq 1) then begin 
    titles[k] = 'Mnfw'
    k=k+1  
endif
; Conc
if ((fit_type[2] eq 1)) then begin
    titles[k] = 'Conc'
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
    titles[k] = 'Roff'
endif

;-------------------------------------------------------------------------
; Make a histogram of marginal posteriors versus the Priors
;-------------------------------------------------------------------------

if keyword_set(hist) then begin

    !p.multi=[0,npars,1]   ; Default   

    if (npars eq 2) then begin
        !p.multi=[0,2,1]        ; scale to number of params
    endif
    if (npars eq 3 or npars eq 4) then begin
        !p.multi=[0,2,2]        ; scale to number of params
    endif

    xx=findgen(100)/99.0
 
    ; Go through the parameters
    for i=0, npars-1 do begin

        bin = p_sigma[i]/5.0
        xr  = long(100*p_mean[i])/100.0+long(10.0*p_sigma[i]*[-4,4])/10.0
        
        plothist,reform(pars[i,*]),bin=bin,hx,hy,/noplot  ; Histogram
        ymax=1000*ceil(max(hy)*1.1/1000.0)

        plothist,reform(pars[i,*]),bin=bin,_extra=ex,hx,hy,xr=xr,xticks=2,yticks=2,xtit=titles[i],/xst,yr=[0,ymax],/yst
        max_p=max(pars[i,*])  ;max of chain
        min_p=min(pars[i,*])  ;min of chain
        xr=(max_p-min_p)
            
        ;plot the prior
        xxx=findgen(1000)/999*8-4                ; [-4,4]
        uuu=(xxx-Prior_mean[i])/Prior_sigma[i]   ; Mean and sigma of prior
        yyy=exp(-0.5*uuu^2)                      ; Gaussian prior
        yyy=(yyy/max(yyy))*max(hy)*1.05          ; Normalize for the plot
        oplot,xxx,yyy,color=!red

    endfor
   
;-------------------------------------------------------------------------
; Marginal Posterior distribution for parameter pairs
;------------------------------------------------------------------------- 

endif else begin
    
    !p.multi=[0,npars,npars]
    !p.thick=3
    !x.thick=3
    !y.thick=3
    !p.charsize=1.2
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

    for ii=0, npars-1 do begin         ; Go through params
        for jj=ii+1, npars-1 do begin     ; Go through params


           blank=replicate(' ',30)
           fill=replicate('',30)
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
                
;                xp=xp+xran
;                yp=yp+yran
                
           sm=7                 ; width of gaussian smoothing
                
                ; NOTE : this was a conflict between hist2d in the ICG dir and dave's version
                ; Now using the version that is in /idl/Alexie

;           ploth, xp, yp, h, psym=3,/yno,_extra=ex,xtit=xtit,ytit=ytit,g_sm=sm,/noplot ; h is the result

;           density_contours,xp,yp,h,lev=lev,/noplot
   
;           l=lindgen(100*100L)
;           xa=float(reform(l mod 100,100,100))
;           ya=float(reform(l / 100,100,100))
;           xr=long(100*p_mean[i])/100.0+long(100.0*p_sigma[i]*[-4,4])/100.0
;           yr=long(100*p_mean[j])/100.0+long(100.0*p_sigma[j]*[-4,4])/100.0
;           contour,h.map,xa*(h.xrange[1]-h.xrange[0])/100.0+h.xrange[0],ya*(h.yrange[1]-h.yrange[0])/100.0+h.yrange[0],xtit=xtit,ytit=ytit,xticks=2,$
;                   yticks=2,levels=lev,/fill,c_color=[!red,!blue],xr=xr,yr=yr,/xst,/yst,position=[x1,y1,x2,y2]
;           print,nn

           rangeFactor=3.
           xr=p_mean[ii]+rangeFactor*p_sigma[ii]*[-1,1]
           yr=p_mean[jj]+rangeFactor*p_sigma[jj]*[-1,1]
           nBins=50.
           xbin=(xr[1]-xr[0])/nBins
           ybin=(yr[1]-yr[0])/nBins
           plot,/nodata,xr,yr,xst=1,yst=1,xtit=xtit,ytit=ytit,position=[x1,y1,x2,y2],xtickn=xtickn,ytickn=ytickn
           oplot_contours,xp,yp,xbin,ybin,xr,yr,!black;,g_smooth=sm
        endfor
    endfor
    
    device,/close

endelse


return
end



