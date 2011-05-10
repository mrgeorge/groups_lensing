; copied from /Users/alexie/idl/MCMC/dsmodel/plot_halofit.pro on July 12 2010
; add NFW offset model

pro plot_halofit_offset, zarray, file, ps_file, fit_t, log_mass, log_mass_err, p_mean, p_sigma, title,$
                         makeplot=makeplot,$
                         quick_c=quick_c,$
                         pbzk=pbzk,$
                         show_result = show_result,$
                         dont_plot=dont_plot,$ 
                         cov=cov,$
                         other_fit=other_fit,$
                         groups=groups,$
                         jackknife=jackknife,$
                         fast=fast,$ ; this controls MCMC nstep and burning
                         panelplot=panelplot,$
                         use_200=use_200,$
                         maccio=maccio,$
                         weak_shear_corr=weak_shear_corr,$
                         bcg=bcg,$
                         xmarge=xmarge,$
                         ymarge=ymarge,$
                         xch=xch,$
                         ych=ych,$
                         notitle=notitle,$
                         agn=agn,$
                         stackx=stackx,$ ; stack the x axis by r/r200
                         center=center,$
                         refcen=refcen

; ---  For the moment just using regular errors
; ---  Can also switch to full cov matrix if necessary. 

start=systime(1)

;-------------------------------------------------------------------------
; Program that calls the MCMC routine to fit data and then makes a plot.
; Note: Mass is virial.
;-------------------------------------------------------------------------

common fit_options, q_c, lens_redshift, fit_type, lens_m_sun, log_sm, use_m200, neg_points, pos_points,str2,str3, use_group, use_maccio, xmar, ymar, xchars, ychars, no_title, ws_corr, lz_mean, sx

;-------------------------------------------------------------------------
; Read Structure
; ** NOTE : NEED TO CHEK THE DIFFERENCE BETWEEN MEAN AND MEDIAN
;-------------------------------------------------------------------------

full_str = mrdfits(file,1)

    ; Testing the ZSPEC correction here:
    ; full_str.we1_mean = full_str.we1_mean*1.07

sel_points     = where(full_str.e1_num ge 15)
sel_points_neg = where(full_str.e1_num ge 15 and full_str.we1_mean le 0)
sel_points_pos = where(full_str.we1_mean gt 0)

if keyword_set(groups) then begin
    ; don't use the innermost bin
    ;sel_points     = where(full_str.e1_num ge 15 and full_str.plot_radius_kpc gt 50)  ; this is for the fit
    ;sel_points_neg = where((full_str.we1_mean le 0 or full_str.e1_num lt 15) and full_str.plot_radius_kpc gt 50)
    ;sel_points_pos = where(full_str.we1_mean gt 0 and full_str.plot_radius_kpc gt 50)

    ; Changed here for Matt's centering paper
    if keyword_set(stackx) then begin
        sel_points     = where(full_str.e1_num ge 10) ; this is for the fit
        sel_points_neg = where(full_str.we1_mean le 0 and full_str.e1_num gt 10)
        sel_points_pos = where(full_str.we1_mean gt 0 and full_str.e1_num gt 10)
    endif else begin
        sel_points     = where(full_str.e1_num ge 10 and full_str.plot_radius_kpc gt 20) ; this is for the fit
        sel_points_neg = where((full_str.we1_mean le 0 or full_str.e1_num lt 10) and full_str.plot_radius_kpc gt 20)
        sel_points_pos = where(full_str.we1_mean gt 0 and full_str.plot_radius_kpc gt 20)
    endelse
endif

str={                                            $
meanr:full_str.plot_radius_kpc(sel_points)      ,$
sigma:full_str.we1_mean(sel_points)             ,$
sigmaerr:full_str.we1_error(sel_points)         ,$
sigmaerr_jack:full_str.j_knife_err(sel_points)  ,$
z_lens: full_str.z_lens                         ,$
mean_lz: full_str.mean_lz                       ,$
mean_lz_test: full_str.mean_lz_test             ,$
msun_lens: full_str.msun_lens          }

neg_points     = 0
if (sel_points_neg[0] ne -1) then begin
    neg_points = 1
    str2={                                                 $
         meanr:full_str.plot_radius_kpc(sel_points_neg)   ,$
         sigma:full_str.we1_mean(sel_points_neg)          ,$
         sigmaerr:full_str.we1_error(sel_points_neg)  }
endif

pos_points     = 0
if (sel_points_pos[0] ne -1) then begin
    pos_points = 1
    str3={                                                 $
         meanr:full_str.plot_radius_kpc(sel_points_pos)   ,$
         sigma:full_str.we1_mean(sel_points_pos)          ,$
         sigmaerr:full_str.we1_error(sel_points_pos)  }
endif

lz_mean = full_str.mean_lz

q_c = 0
if keyword_set(quick_c) then q_c = 1

use_group = 0
if keyword_set(groups) then use_group = 1

use_m200 = 0
if keyword_set(use_200) then use_m200 = 1

ws_corr = 0
if keyword_set(weak_shear_corr) then ws_corr =1

use_maccio = 0
if keyword_set(maccio) then use_maccio = 1

if keyword_set(xmarge) then xmar = xmarge else xmar=[-99,-99]

if keyword_set(ymarge) then ymar = ymarge else ymar=[-99,-99]

if keyword_set(xch) then xchars=xch else xchars=1.0

if keyword_set(ych) then ychars=ych else ychars=1.0

if keyword_set(notitle) then no_title=notitle else no_title=0

sx=0
if keyword_set(stackx) then sx=1

lens_redshift = full_str.z_lens

fit_type = fit_t

; point mass
log_sm          = str.msun_lens
lens_m_sun      = (10^(log_sm))/1e12   ; 10^12 h^1 Msun

;-------------------------------------------------------------------------
; PLOT STUFF
;-------------------------------------------------------------------------

;yr = [0.1,1000]
yr = [0.3,1000]
xr = [0.01,4]

if keyword_set(pbzk) then xr = [0.005,6] 
if keyword_set(groups) then xr = [0.02,2];xr = [0.02,6] ;xr = [0.04,6] ;[0.03,6] 
if keyword_set(groups) then yr = [0.5,3000] 
if keyword_set(agn) then xr=[0.01,2.5]

;-------------------------------------------------------------------------
; Number of Params to fit
;-------------------------------------------------------------------------

select = where(fit_type eq 1)
npars = n_elements(select)

;-------------------------------------------------------------------------
; MCMC METHOD
;-------------------------------------------------------------------------

if keyword_set(fast) then begin  ;for quick debugging
    nstep=4000
    burnin=100
endif else begin
    nstep=50000L
    burnin=500L
endelse
    
print,nstep
print,burnin

nch = nstep-burnin  ; Number of elements in the chain

; Call the MCMC procedure. NOTE : 'pars' is the MCMC chain
ds_mcmc, str, pars, p_mean=p_mean, p_sigma=p_sigma, rob_p_mean=rob_p_mean, rob_p_sigma=rob_p_sigma, nstep=nstep, burnin=burnin, center=center,refcen=refcen

; Robust estimation versus regular estimation
print,'----------'
print,'COMPARE regular versus rob :'
print,p_mean,rob_p_mean
print,p_sigma,rob_p_sigma
print,'----------'

print,'END OF DS_MCMC ROUTINE'

if keyword_set(show_result) then begin
      print_result,p_mean,str.z_lens
  endif
print,'-------'
;print,p_mean
;print,p_sigma
;print,rob_p_mean
;print,rob_p_sigma
; stop        
;p_mean     = [13.589287]
;p_sigma    = [0.10999367]
;rob_p_mean = [13.589315]
;rob_p_sigma= [0.10994238]

;-------------------------------------------------------------------------
; MAKE PLOTS
;-------------------------------------------------------------------------

; >>> Data plot >>>
if keyword_set(makeplot) then begin
ps_open2,/PORTRAIT,/ENCAPSULATED,XSIZE=7.0,YSIZE=5,/COLOR,THICK=3,/PS_FONT
device, /cmyk, /times,font_size=14  ;/bold
endif

    if keyword_set(groups) then begin
    make_plot,zarray,str,full_str, p_mean, p_sigma, title,xr=xr,yr=yr,tit=tit,xnoedge=xnoedge,ynoedge=ynoedge,mcmc=mcmc,/groups,center=center,refcen=refcen
    endif else begin
    make_plot,zarray,str,full_str, p_mean, p_sigma, title,xr=xr,yr=yr,tit=tit,xnoedge=xnoedge,ynoedge=ynoedge,mcmc=mcmc
    endelse


if keyword_set(makeplot) then begin
ps_close
file1 = strcompress(ps_file+'.eps',/remove_all)
command='mv idl.eps '+file1
spawn,command
endif


if (NOT keyword_set(panelplot)) then begin

; >>> Cov plot1 >>>
if keyword_set(makeplot) then begin
ps_open2,/PORTRAIT,/ENCAPSULATED,XSIZE=7.0,YSIZE=5,/COLOR,THICK=3,/PS_FONT
device, /cmyk, /times,font_size=14  ;/bold
endif

  ds_cov_plots, pars, p_mean=p_mean, p_sigma=p_sigma  

if keyword_set(makeplot) then begin
ps_close
file1 = strcompress(ps_file+'_cov.eps',/remove_all)
command='mv idl.eps '+file1
spawn,command
endif

; >>>> Cov plot2 (with histogram) >>>
;if keyword_set(makeplot) then begin
;ps_open2,/PORTRAIT,/ENCAPSULATED,XSIZE=10.0,YSIZE=5,/COLOR,THICK=3,/PS_FONT
;device, /cmyk, /times,font_size=14  ;/bold
;endif

;  ds_cov_plots, pars, p_mean=p_mean, p_sigma=p_sigma,/hist  

;if keyword_set(makeplot) then begin
;ps_close
;file1 = strcompress(ps_file+'_cov_hist.eps',/remove_all)
;command='mv idl.eps '+file1
;spawn,command
;endif

; >>> Convergence plot (convergence of MCMC chain)

;if keyword_set(makeplot) then begin
;ps_open2,/PORTRAIT,/ENCAPSULATED,XSIZE=10.0,YSIZE=7,/COLOR,THICK=3,/PS_FONT
;device, /cmyk, /times,font_size=14  ;/bold
;endif

;group_mcmc_check_conv, pars, p_mean, p_sigma, rob_p_mean, rob_p_sigma

;if keyword_set(makeplot) then begin
;ps_close
;file1 = strcompress(ps_file+'_conver.eps',/remove_all)
;command='mv idl.eps '+file1
;spawn,command
;endif

endif ; panelplot

;-------------------------------------------------------------------------
; MASS
; !!!! Use the robust estimation here !!!!
; (in noisy cases, the mean can be different from the robust mean)
;-------------------------------------------------------------------------

; Baryonic term
; Error on nfw ( a = exp(b) -> da = a.db )
if (fit_type[0] eq 1) then begin   ; Fit for the baryons (r is now p_mean[1])
    log_mass     = rob_p_mean[1]
    log_mass_err = rob_p_sigma[1]
endif else begin
    log_mass     = rob_p_mean[0]
    log_mass_err = rob_p_sigma[0]
endelse

;-------------------------------------------------------------------------
; Make the Correlation plot
;-------------------------------------------------------------------------

;if not keyword_set(groups) then begin
;corr_file = strcompress(file+'.corr.sav',/REMOVE_ALL)

;restore,corr_file

;ps_open2,/PORTRAIT,/ENCAPSULATED,XSIZE=7.0,YSIZE=5,/COLOR,THICK=3,/PS_FONT
;tvim2,corr,/scale
;ps_close

;endif

;command='mv idl.eps '+ps
;spawn,command

;-------------------------------------------------------------------------
; Convergence
;-------------------------------------------------------------------------

print,''
print,'Runtime = ',(systime(1)-start)/3600.0, ' Hours' 
print,'Runtime = ',(systime(1)-start)/60.0, ' Minutes' 

return
end
