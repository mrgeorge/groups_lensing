pro run_ds_mcmc, lens_infile, $
                 fit_t, $
                 rob_p_mean, rob_p_sigma, $
                 fast=fast, medium=medium, slow=slow, $
                 stackx=stackx, $
                 chainFile=chainFile,$
                 burnin=burnin, $
                 noSave=noSave, $
                 ps=ps, sis=sis, tis=tis, rhotis=rhotis, $
                 off2dDelta=off2dDelta, off3dDelta=off3dDelta, off3dMax=off3dMax

; Partial replacement of plot_halofit - only does the modeling, plotting is done elsewhere

;-------------------------------------------------------------------------
; Program that calls the MCMC routine to fit data and then makes a plot.
; Note: Mass is virial.
;-------------------------------------------------------------------------

common fit_options, q_c, lens_redshift, fit_type, lens_m_sun, log_sm, use_m200, neg_points, pos_points,str2,str3, use_group, use_maccio, xmar, ymar, xchars, ychars, no_title, ws_corr, lz_mean, sx, cen_type, off_type


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

; FILL COMMON BLOCK VARIABLES FOR DS_MCMC
fit_type = fit_t
lz_mean = full_str.mean_lz
lens_redshift = full_str.z_lens
log_sm = str.msun_lens
lens_m_sun = (10^(log_sm))/1.e12 ; 10^12 h^1 Msun

sx = keyword_set(stackx)        ; /stackx

; model for central component
if(fit_t[0] NE 0) then begin
   if(keyword_set(ps)) then $
      cen_type='ps' $
   else if(keyword_set(sis)) then $
      cen_type='sis' $
   else if(keyword_set(tis)) then $
      cen_type='tis' $
   else if(keyword_set(rhotis)) then $
      cen_type='rhotis' $
   else begin
      print,'RUN_DS_MCMC: must choose a cen_type (ps|sis|tis|rhotis)'
      stop
   endelse
endif

; offset model
if(fit_t[6] NE 0) then begin
   if(keyword_set(off2dDelta)) then $
      off_type='delta2d' $
   else if(keyword_set(off3dDelta)) then $
      off_type='delta3d' $
   else if(keyword_set(off3dMax)) then $
      off_type='max3d' $
   else begin
      print,'RUN_DS_MCMC: must choose an off_type (off2dDelta|off3dDelta|off3dMax)'
      stop
   endelse
endif

; hardcoding options here for other program calls
; this is for stacking xray groups and studying centers
q_c = 1                         ; /quick_c
use_group = 1                   ; /groups
use_m200 = 1                    ; /use_200
ws_corr = 0                     ; /weak_shear_corr
use_maccio = 0                  ; /use_maccio


;-------------------------------------------------------------------------
; MCMC METHOD
;-------------------------------------------------------------------------

if keyword_set(fast) then begin  ;for quick debugging
    nstep=4000
    burnin=100
endif else if(n_elements(medium) GT 0) then begin
    nstep=15000L
    burnin=200L
endif else if(n_elements(slow) GT 0) then begin
    nstep=1000000L
    burnin=1000L
endif else begin
    nstep=50000L
    burnin=500L
endelse
    
print,'MCMC Nstep ',nstep
print,'MCMC burn-in ',burnin

nch = nstep-burnin  ; Number of elements in the chain

; Call the MCMC procedure. NOTE : 'pars' is the MCMC chain
ds_mcmc, str, pars, p_mean=p_mean, p_sigma=p_sigma, rob_p_mean=rob_p_mean, rob_p_sigma=rob_p_sigma, nstep=nstep, burnin=burnin, chainFile=chainFile

; Robust estimation versus regular estimation
print,'----------'
print,'COMPARE regular versus rob :'
print,p_mean,rob_p_mean
print,p_sigma,rob_p_sigma
print,'----------'

if(n_elements(noSave) EQ 0) then begin
   ; Write fit type and parameters to struct/file
   if(NOT(tag_exist(full_str,'FIT_TYPE'))) then $
      full_str=create_struct(full_str,'FIT_TYPE',fit_type,'P_MEAN',rob_p_mean,'P_SIGMA',rob_p_sigma) $
   else if(NOT(tag_exist(full_str,'FIT_TYPE2'))) then $
      full_str=create_struct(full_str,'FIT_TYPE2',fit_type,'P_MEAN2',rob_p_mean,'P_SIGMA2',rob_p_sigma) $
   else begin
      print,'run_ds_mcmc: tags FIT_TYPE and FIT_TYPE2 already exist in this struct' ; it's already been modeled twice, if you want more than two models rethink the data structure and code
      stop
   endelse

   mwrfits,full_str,lens_infile,/create
endif

print,'END OF DS_MCMC ROUTINE'

; rob_p_mean and rob_p_sigma are filled and returned and written to lens_infile
end
