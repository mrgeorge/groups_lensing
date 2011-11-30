pro ds_cov_plots, pars, p_mean=p_mean, p_sigma=p_sigma, hist=hist

; copied from /Users/alexie/idl/MCMC/dsmodel/ds_cov_plots.pro with changes

;-----------------------------------------------------------------------------
; Adapted from Dave : /Users/alexie/idl/dave_idl/inversion/mcmc_cov_plots3.pro
;-----------------------------------------------------------------------------

common mcmc_block, x, y, ivar, psigma, npars, Prior_mean, Prior_sigma
common fit_options, q_c, lens_redshift, fit_type, lens_m_sun

;-- Get robust measures
mcmc_stats,pars,p_mean,p_sigma,p_mean_rob,p_sigma_rob,nsig=nsig,niter=niter

sz=size(pars,/dim)   ; pars = chain
npars=sz[0]              ; n=number of pars

;!p.charsize=2

;!x.omargin=[4,2]
;!y.omargin=[3,2]

;!x.margin=[7,3]
;!y.margin=[3,2]

; -- Set up plot Labels -------------------------------------------------
; 0  baryonic mass
; 1  R_vir : NFW virial radius
; 2  C     : NFW concentration
; 3  alpha : fraction
; 4  Bias
; 5  m_sigma : dispersion in central mass
;-------------------------------------------------------------------------

p_lab = strarr(npars)

k = 0

; Baryonic term
if (fit_type[0] eq 1) then begin
    p_lab[k] = 'Mcen'
    k=k+1  
endif
; Rvir
if (fit_type[1] eq 1) then begin 
    p_lab[k] = 'Mnfw'
    k=k+1  
endif
; Conc
if ((fit_type[2] eq 1)) then begin
    p_lab[k] = 'Conc'
    k=k+1  
endif
; Alpha
if (fit_type[3] eq 1) then begin 
    p_lab[k] = 'Alpha'
    k=k+1  
endif
; Bias
if (fit_type[4] eq 1) then begin 
    p_lab[k] = 'Bias'
    k=k+1  
endif
; Dispersion in central mass
if (fit_type[5] eq 1) then begin 
    p_lab[k] = 'Dispersion'
    k=k+1  
endif
; Offset radius
if (fit_type[5] eq 1) then begin 
    p_lab[k] = 'Roff'
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

        bin = p_sigma_rob[i]/5.0
        xr  = long(100*p_mean[i])/100.0+long(10.0*p_sigma[i]*[-4,4])/10.0
        
        plothist,reform(pars[i,*]),bin=bin,hx,hy,/noplot  ; Histogram
        ymax=1000*ceil(max(hy)*1.1/1000.0)

        plothist,reform(pars[i,*]),bin=bin,_extra=ex,hx,hy,xr=xr,xticks=2,yticks=2,xtit=p_lab[i],/xst,yr=[0,ymax],/yst
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
    
    !p.multi=[0,1,0]
    
    for i=0, npars-1 do begin         ; Go through params
        for j=0, npars-1 do begin     ; Go through params

            if i lt j then begin

                xtit=p_lab[i]  ;title
                ytit=p_lab[j]  ;title

                xw=sigma(pars[i,*])
                yw=sigma(pars[j,*])
                nn=n_elements(pars[i,*])
                xp=pars[i,*]
                yp=pars[j,*]
                xran=randomn(seed,nn)*xw/5.0
                yran=randomn(seed,nn)*yw/5.0
                
;                xp=xp+xran
;                yp=yp+yran
                
                sm=7  ; width of gaussian smoothing
                
                ; NOTE : this was a conflict between hist2d in the ICG dir and dave's version
                ; Now using the version that is in /idl/Alexie

                ploth, xp, yp, h, psym=3,/yno,_extra=ex,xtit=xtit,ytit=ytit,g_sm=sm,/noplot    ; h is the result

                density_contours,xp,yp,h,lev=lev,/noplot
   
                l=lindgen(100*100L)
                xa=float(reform(l mod 100,100,100))
                ya=float(reform(l / 100,100,100))
                xr=long(100*p_mean[i])/100.0+long(100.0*p_sigma[i]*[-4,4])/100.0
                yr=long(100*p_mean[j])/100.0+long(100.0*p_sigma[j]*[-4,4])/100.0
                contour,h.map,xa*(h.xrange[1]-h.xrange[0])/100.0+h.xrange[0],ya*(h.yrange[1]-h.yrange[0])/100.0+h.yrange[0],xtit=xtit,ytit=ytit,xticks=2,$
                  yticks=2,levels=lev,/fill,c_color=[!red,!blue],xr=xr,yr=yr,/xst,/yst
                print,nn
            endif
        endfor
    endfor
    
endelse


return
end



