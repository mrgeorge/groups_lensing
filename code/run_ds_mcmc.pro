pro run_ds_mcmc, lens_infile, $
                 fit_t, $
                 rob_p_mean, $
                 rob_p_sigma, $
                 fast=fast

; Partial replacement of plot_halofit - only does the modeling, plotting is done elsewhere

;-------------------------------------------------------------------------
; Program that calls the MCMC routine to fit data and then makes a plot.
; Note: Mass is virial.
;-------------------------------------------------------------------------

common fit_options, q_c, lens_redshift, fit_type, lens_m_sun, log_sm, use_m200, neg_points, pos_points,str2,str3, use_group, use_maccio, xmar, ymar, xchars, ychars, no_title, ws_corr, lz_mean, sx


; Read structure for measured lensing signal
full_str = mrdfits(lens_infile,1)

; Only include bins with >=10 background sources
sel_points = where(full_str.e1_num GE 10)

str={                                            $
meanr:full_str.plot_radius_kpc[sel_points]      ,$
sigma:full_str.we1_mean[sel_points]             ,$
sigmaerr:full_str.we1_error[sel_points]         ,$
sigmaerr_jack:full_str.j_knife_err[sel_points]  ,$
z_lens: full_str.z_lens                         ,$
mean_lz: full_str.mean_lz                       ,$
mean_lz_test: full_str.mean_lz_test             ,$
msun_lens: full_str.msun_lens          }

; Fill some common block variables
fit_type = fit_t
lz_mean = full_str.mean_lz
lens_redshift = full_str.z_lens
log_sm = str.msun_lens
lens_m_sun = (10^(log_sm))/1.e12 ; 10^12 h^1 Msun

; hardcoding options here for other program calls
; this is for stacking xray groups and studying centers
q_c = 1                         ; /quick_c
use_group = 1                   ; /groups
use_m200 = 1                    ; /use_200
ws_corr = 0                     ; /weak_shear_corr
use_maccio = 0                  ; /use_maccio
sx = 1                          ; /stackx

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
    
print,'MCMC Nstep ',nstep
print,'MCMC burn-in ',burnin

nch = nstep-burnin  ; Number of elements in the chain

; Call the MCMC procedure. NOTE : 'pars' is the MCMC chain
ds_mcmc, str, pars, p_mean=p_mean, p_sigma=p_sigma, rob_p_mean=rob_p_mean, rob_p_sigma=rob_p_sigma, nstep=nstep, burnin=burnin

; Robust estimation versus regular estimation
print,'----------'
print,'COMPARE regular versus rob :'
print,p_mean,rob_p_mean
print,p_sigma,rob_p_sigma
print,'----------'

print,'END OF DS_MCMC ROUTINE'

; rob_p_mean and rob_p_sigma are filled and returned
end
