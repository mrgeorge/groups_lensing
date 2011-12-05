pro mcmc_setup, s

;-------------------------------------------------------------------------
; *** SETUP ****
; Set up the priors
;-------------------------------------------------------------------------

common mcmc_block, x, y, ivar, psigma, npars, Prior_mean, Prior_sigma
common fit_options, q_c, lens_redshift, fit_type, lens_m_sun, log_sm, use_m200, neg_points, pos_points,str2,str3, use_group, use_maccio, xmar, ymar, xchars, ychars, no_title, ws_corr, lz_mean, sx,sis,tsis

;-------------------------------------------
; Prior Gaussian MEANS and Gaussian Sigma
; uses logarithmic variables for most except pc, pc= 1/(1+exp(-q))
;-------------------------------------------

if (use_group eq 1) then begin ; Groups
    Prior_mean = setup_4(/print_info)
    Prior_sigma = setup_4(/sigma)
    npars = n_elements(Prior_mean)
    print,'** Number of free parameters : ', npars
    print,'GROUPS'
endif else begin     ; Galaxies
; NEED TO CHECK THE PRIOIRS HERE
    Prior_mean = setup_2(/print_info)
    Prior_sigma = setup_2(/sigma)
    npars = n_elements(Prior_mean)
    print,'** Number of free parameters : ', npars
endelse

;-------------------------------------------
; Struct
; Check error units here
;-------------------------------------------

x     = s.meanr/1e3
y     = s.sigma
yerr  = s.sigmaerr
ivar  = 1.0/yerr^2

; -> Set this up to use the full covariance matrix
;if 1 then begin
;    cov=s.covariance
;    fix_cov,cov,covf
;    ivar=invert(covf)
;endif else begin
;    ivar=1.0/yerr^2
;endelse

;-------------------------------------------
; Run ds_model once just to set up xi_lin
;-------------------------------------------
; MRG - this line commented in order to switch to using get_ds_model
;       ONLY, replacing ds_model. This will cause modeling of the
;       two-halo term to fail, since it's currently not
;       supported by get_ds_model.
;fake_model = ds_model(x,Prior_mean)

;-------------------------------------------
; Sample from the prior
; psigma is just the Guassian distribution used to make the step size. 
; This should be smaller than the prior width so I shrink it. 
; The efficiency (speed of convergence) of MCMC depends on getting the
; step size not too big and not too small. 
;-------------------------------------------

; This should have the same number of elements as nparams
shrink = fltarr(n_elements(Prior_sigma))
shrink(*) = 0.5

psigma = Prior_sigma*shrink

return
end

function mcmc_check_limits, res
; return 0 if pars are within bounds, 1 if it hits any limits

common fit_options, q_c, lens_redshift, fit_type, lens_m_sun, log_sm, use_m200, neg_points, pos_points,str2,str3, use_group, use_maccio, xmar, ymar, xchars, ychars, no_title, ws_corr, lz_mean, sx,sis,tsis

   ii=0
   if(fit_type[0] EQ 1) then begin ; Mcen - require 10<Mcen<13
      if(res[ii] LE 10 OR res[ii] GE 13) then return, 1
      ii+=1
   endif
   if(fit_type[1] EQ 1) then begin ; Halo mass - require 10<Mhalo<16
      if(use_group EQ 1 AND res[ii] LE 10.) then return, 1 $
      else if(use_group EQ 1 AND res[ii] GE 16.) then return, 1
      ii+=1
   endif
   if(fit_type[2] EQ 1) then begin ; concentration - require 1<conc<10
      if(res[ii] LE 1. OR res[ii] GE 10.) then return, 1
      ii+=1
   endif
   if(fit_type[3] EQ 1) then begin ; alpha
      ii+=1
   endif
   if(fit_type[4] EQ 1) then begin ; bias
      ii+=1
   endif
   if(fit_type[5] EQ 1) then begin ; sigma_Mh
      ii+=1
   endif
   if(fit_type[6] EQ 1) then begin ; offset - require 0<Roff<1Mpc
      if(res[ii] LE 0. OR res[ii] GE 1.) then return, 1
      ii+=1
   endif

return, 0
end

;-------------------------------------------------------------------------
; *** STEP ***
;-------------------------------------------------------------------------

function mcmc_step, seed, pars

common mcmc_block, x, y, ivar, psigma, npars, Prior_mean, Prior_sigma
common fit_options, q_c, lens_redshift, fit_type, lens_m_sun, log_sm, use_m200, neg_points, pos_points,str2,str3, use_group, use_maccio, xmar, ymar, xchars, ychars, no_title, ws_corr, lz_mean, sx,sis,tsis
   ;Gaussian steps, two choices
   ;step away from previous set of pars

   repeat begin
      res=psigma*randomn(seed,npars)+pars ; ! Note : this is where psigma comes into play !
   endrep until (mcmc_check_limits(res) EQ 0)
   
return,res
end 


;-------------------------------------------------------------------------
; *** LIKELIHOOD ***
;-------------------------------------------------------------------------

function mcmc_like, p

common mcmc_block, x, y, ivar, psigma, npars, Prior_mean, Prior_sigma
common fit_options, q_c, lens_redshift, fit_type, lens_m_sun, log_sm, use_m200, neg_points, pos_points,str2,str3, use_group, use_maccio, xmar, ymar, xchars, ychars, no_title, ws_corr, lz_mean, sx,sis,tsis

  ; reject point if it is out of bounds by returning -Inf likelihood
  if(mcmc_check_limits(p)) then return, -!values.f_infinity


  ; Inverse variance (ivar=1.0/yerr^2)
  ; Basically, this is the number of data points
  n_dim = size(ivar, /n_dim)

  ;---------------------------------------------------
  ; THIS IS THE DELTA SIGMA MODEL : model=delta sigma
  ;--------------------------------------------------

;  model = ds_model(x,p,/skip_dslin)
  get_ds_model,fit_type,p,lens_redshift,log_sm,x,tot=model,sis=sis,tsis=tsis,use_m200=use_m200
  
  ; check dave's version which is : nfw_delta_sigma_mcmc(x,p)

  diff = y-model  ;(array of X_i-mu_i)

  ; This computes the chi square = Sum (X_i-mu_i)^2/sigma_i^2
  if n_dim eq 1 then begin 
      chi2 = total(diff^2*ivar )
  endif else begin 
      chi2 = transpose(diff)#(reform(ivar##diff))
  endelse 

  chi_square = chi2
  
  ;if 1 then begin
  ;    ; This also include the chi square with respect to initial guess
  ;    chi2_prior=total(((p-Prior_mean)/Prior_sigma)^2)
  ;    chi_square=chi2+chi2_prior
  ;endif
  
  ; This is the log likelihood (-1/2 * (sum(x-mod)/sigma)^2 )
  return, -0.5*chi_square
  
end 

;-------------------------------------------------------------------------
; *** MAIN PROGRAM ***
; p_sigma      = Best fit parameters from chain
; p_mean       = Best fit parameters from chain
; npars        = number of free parameters
; pars         = Markov chain : [n_params, nsteps]
; ivar         = inverse variance of signal (ivar=1.0/yerr^2)
; parguess     = parameter guess
; Prior_mean   = what you put in as a prior
; Prioir_sigma = what you put in as sigma
;-------------------------------------------------------------------------

pro ds_mcmc,s, pars, p_mean=p_mean, p_sigma=p_sigma, rob_p_mean=p_mean_rob,$
            rob_p_sigma=p_sigma_rob,$
            lum=lum,$
            zoom=zoom,$
            nstep=nstep,$
            burnin=burnin,$
            rand_start=rand_start,$
            chainFile=chainFile

if n_params() eq 0 then begin
    print,'-syntax mcmc_nfw_conc,s,pars,p_mean=p_mean,p_sigma=p_sigma,plot=plot,lum=lum'
    return
endif

common mcmc_block, x, y, ivar, psigma, npars, Prior_mean, Prior_sigma
common fit_options, q_c, lens_redshift, fit_type, lens_m_sun, log_sm, use_m200, neg_points, pos_points,str2,str3, use_group, use_maccio, xmar, ymar, xchars, ychars, no_title, ws_corr, lz_mean, sx,sis,tsis

; Number of steps
if n_elements(nstep) eq 0 then nstep = 10000
nstep_quick=min([long(nstep/3.0),10000L])  ; this is for first run
if n_elements(burnin) eq 0 then burnin = 500

if nstep lt burnin*1.5 then begin
    print,'ERROR - nstep should be much bigger than burnin'
    stop
endif

mcmc=obj_new('mcmc')

;--- Call Setup Function ---
;--- Defines the Priors (Prior_mean and Prior_sigma)
mcmc_setup,s

;--- Runs twice.....
;--- Set the initial guess to the prioir
parguess = Prior_mean

;--- Start Randomly
if keyword_set(rand_start) then begin
    print,'Random start'
    parguess=parguess+randomn(seed,npars)*psigma*3
endif

;--- The first time is with nstep_quick
print,'>> First MCMC run'
pars = mcmc->run('mcmc_step', 'mcmc_like', $
                 nstep_quick, parguess, printstep=10000,/log)

;-- ????
count_trans,pars,new,k
print,'Num transitions ',k
trans=nstep_quick/float(k)
print,'Average Transition length ', trans, '  steps'

;--- get rid of burnin part before calculating mean and sigma
pars=pars[*,burnin:*]

; --- Calculate mean and sigma with clipping
; --- This is where robust mean and sigma are calculated (p_mean_rob, p_sigma_rob)
; --- nsig=sigma for sigma clipping, niter=for sigma clipping
mcmc_stats,pars,p_mean,p_sigma,p_mean_rob,p_sigma_rob,nsig=4,niter=5

; Why compare to psigma here ..??????
prat=p_sigma_rob/psigma
print
print,'Psig ratio1',prat
print

;-- Now use starting points from previous run --
parguess = p_mean_rob
psigma   = p_sigma_rob

;-- Second Run
print,'>> Second MCMC run'
pars = mcmc->run('mcmc_step', 'mcmc_like', $
                 nstep, parguess, printstep=10000,/log,file=chainFile)
if(n_elements(chainFile) GT 0) then pars=mcmc->read_trials(chainFile)

; ??? what is k ??
count_trans,pars,new,k
print,k,' transitions'
kmin=10
if k lt kmin then begin
    print,'ERROR - number of transitions too small',k
    print,'Change psigma'
    stop
endif

; ???? should this be nstep instead ????
trans=nstep_quick/float(k)
print,'Average Transition length ', trans, '  steps'

;--- get rid of burning part before calculating mean and sigma
pars=pars[*,burnin:*]

; --- Calculate mean and sigma with clipping
; --- This is where robust mean and sigma are calculated (p_mean_rob, p_sigma_rob)
; --- nsig=sigma for sigma clipping, niter=for sigma clipping
mcmc_stats,pars,p_mean,p_sigma,p_mean_rob,p_sigma_rob,nsig=4,niter=5

prat=p_sigma_rob/psigma
print
print,'Psig ratio2',prat
print

return
end
















