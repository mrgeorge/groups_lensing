function setup_4, sigma=sigma, print_info=print_info

; copied from /Users/alexie/idl/MCMC/dsmodel/setup_4.pro with changes
; to priors

;-------------------------------------------------------------------------
; Do the fit in MASS instead of R200
;-------------------------------------------------------------------------

common fit_options, q_c, lens_redshift, fit_type, lens_m_sun, log_sm, use_m200, neg_points, str2

;-------------------------------------------------------------------------
; MEAN
;-------------------------------------------------------------------------

; 0  M0       : baryonic mass           : !! alog10 !!
; 1  Mnfw     : NFW mass                : !! alog10 !!
; 2  C        : NFW concentration       : linear
; 3  alpha    : fraction of satellites  : q = alog(a/(1-a))
; 4  b        : bias                    : not log
; 5  m_sigma  : dispersion in mass      : See DAVE's paper
; 6  offset   : offset radius           : Mpc

; -- M0
M0      = log_sm

; -- Bias
bias    =  1.0           ; bias

; -- Dispersion in mass
m_sigma =  0.3           ; dispersion in mass selection (log-normal)

; -- Mnfw: ~ 13.5
mnfw = 13.5

; -- Estimated concentration
conc = get_conc(mnfw,lens_redshift,/use200)

; estimated offset
offset = 0.050 ; Mpc

;-------------------------------------------------------------------------
; No Alpha for groups
;-------------------------------------------------------------------------

alpha = 0.1

;-------------------------------------------------------------------------
; SIGMA: Prior Gaussian sigmas
;-------------------------------------------------------------------------

sigma_M0      = 1.0   ; Small : not too constrained
sigma_Mnfw    = 0.8   ; Covers 12.7 -> 14.3 at 1 sigma
sigma_conc    = 3.0   ; ~1 -> 7
sigma_q       = 0.7 
sigma_bias    = 0.3
sigma_m_sigma = 0.1
sigma_offset  = 0.05 ; Mpc

;-------------------------------------------------------------------------
; RETURN PRIORS
;-------------------------------------------------------------------------

beta   =  alpha/(1-alpha) 
q      =  alog(beta)       ;q varies between -3 and 3

Prior_mean = [  $
             M0             ,$  ; alog10 (alreay in this format)
             mnfw           ,$  ; alog10 (alreay in this format)
             conc           ,$  ; conc
             q              ,$  ; This is already in log
             bias           ,$  ; Not in log
             m_sigma        ,$  ; ???
             offset          $
 ]

Prior_sigma = [$
              sigma_M0      ,$
              sigma_Mnfw    ,$
              sigma_conc    ,$
              sigma_q       ,$
              sigma_bias    ,$
              sigma_m_sigma ,$
              sigma_offset   $
]

;-------------------------------------------------------------------------
; PRINT INFO
;-------------------------------------------------------------------------

if (keyword_set(print_info)) then begin
    print,''
    print,'--------------------------------'
    print,'--- SET UP INFORMATION ---------'
    print,''
    print,'Stellar Mass', log_sm
    print,''
    print,'-- NFW TERM --'
    print,'M nfw', mnfw
    print,'Conc', conc

    print,'M min', mnfw-sigma_Mnfw, '    M max ',mnfw+sigma_Mnfw
    print,'Conc min', conc-sigma_conc,'Conc max', conc+sigma_conc
endif

;-------------------------------------------------------------------------
; WHICH PARMETERS ARE ACTUALLY FREE ?
;-------------------------------------------------------------------------

select = where(fit_type eq 1)

Prior_mean  = Prior_mean(select)
Prior_sigma = Prior_sigma(select)

;-------------------------------------------------------------------------
; RETURN PRIORS
;-------------------------------------------------------------------------

if keyword_set(sigma) then begin
    return, Prior_sigma
endif else begin
    return, Prior_mean
endelse

end
