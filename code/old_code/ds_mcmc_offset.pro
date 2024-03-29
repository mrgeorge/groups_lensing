pro mcmc_setup, s

;-------------------------------------------------------------------------
; *** SETUP ****
; Set up the priors
;-------------------------------------------------------------------------

common mcmc_block, x, y, ivar, psigma, npars, Prior_mean, Prior_sigma
common fit_options, q_c, lens_redshift, fit_type, lens_m_sun, log_sm, use_m200, neg_points, pos_points,str2,str3, use_group, use_maccio, xmar, ymar, xchars, ychars, no_title, ws_corr, lz_mean, sx

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

fake_model = ds_model(x,Prior_mean)

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

;-------------------------------------------------------------------------
; *** STEP ***
;-------------------------------------------------------------------------

function mcmc_step, seed, pars

common mcmc_block, x, y, ivar, psigma, npars, Prior_mean, Prior_sigma
common fit_options, q_c, lens_redshift, fit_type, lens_m_sun, log_sm, use_m200, neg_points, pos_points,str2,str3, use_group, use_maccio, xmar, ymar, xchars, ychars, no_title, ws_corr, lz_mean, sx
   ;Gaussian steps, two choices
   ;step away from previous set of pars
   res=psigma*randomn(seed,npars)+pars   ; ! Note : this is where psigma comes into play !

   ; Enforce a reasonable value for the mass for the noisy groups:
   ; (this is imp: prevents runaway results)

   if (fit_type[0] eq 1) then mass=res[1] else mass=res[0]

   ;if (use_group eq 1 AND mass le 11.0) then begin
   if (use_group eq 1 and mass le 10.0) then begin
       if (fit_type[0] eq 1) then res[1]=10.0 else res[0]=10.0
   endif
   if (use_group eq 1 and mass ge 16.0) then begin
       if (fit_type[0] eq 1) then res[1]=16.0 else res[0]=16.0
   endif

   ; NOTE :
   ; This is where you can put a prior on a variable to be more than 0 for example !!!
   ; e.g : NOTE : I MAY NEED TO CHECK THAT m_sigma remains positive
   
return,res
end 

;-------------------------------------------------------------------------
; *** LIKELIHOOD ***
;-------------------------------------------------------------------------

function mcmc_like, p

common mcmc_block, x, y, ivar, psigma, npars, Prior_mean, Prior_sigma
common fit_options, q_c, lens_redshift, fit_type, lens_m_sun, log_sm, use_m200, neg_points, pos_points,str2,str3, use_group, use_maccio, xmar, ymar, xchars, ychars, no_title, ws_corr, lz_mean, sx

  ; Inverse variance (ivar=1.0/yerr^2)
  ; Basically, this is the number of data points
  n_dim = size(ivar, /n_dim)

  ;---------------------------------------------------
  ; THIS IS THE DELTA SIGMA MODEL : model=delta sigma
  ;--------------------------------------------------

  model = ds_model(x,p,/skip_dslin)
  
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

pro ds_mcmc_offset,s, pars, p_mean=p_mean, p_sigma=p_sigma, rob_p_mean=p_mean_rob,$
            rob_p_sigma=p_sigma_rob,$
            lum=lum,$
            zoom=zoom,$
            nstep=nstep,$
            burnin=burnin,$
            rand_start=rand_start,$
            center=center,$
            refcen=refcen

if n_params() eq 0 then begin
    print,'-syntax mcmc_nfw_conc,s,pars,p_mean=p_mean,p_sigma=p_sigma,plot=plot,lum=lum'
    return
endif

common mcmc_block, x, y, ivar, psigma, npars, Prior_mean, Prior_sigma
common fit_options, q_c, lens_redshift, fit_type, lens_m_sun, log_sm, use_m200, neg_points, pos_points,str2,str3, use_group, use_maccio, xmar, ymar, xchars, ychars, no_title, ws_corr, lz_mean, sx
common xGroupCenter, xCenter, xRefCen

if(keyword_set(center)) then xCenter=center else xCenter=''
if(keyword_set(refcen)) then xRefCen=refcen else xRefCen=''

; Number of steps
if n_elements(nstep) eq 0 then nstep = 10000
nstep_quick=long(nstep/3.0)  ; this is for first run
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
                 nstep, parguess, printstep=10000,/log)

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
















