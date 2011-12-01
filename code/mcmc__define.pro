;+
; NAME:
;  MCMC (IDL Class File)
;
; PURPOSE:
;  Calculate a Monte Carlo Markov Chain. 
;
; METHODS:
;  For detailed info for each method, use
;  IDL> doc_method, 'mcmc:methodname'
;
;  mcmc::run(step_func, like_func, nstep, parguess, like=, /log, seed=,
;            printstep=, file=)
;  mcmc::run, step_func, like_func, nstep, parguess, pars, like, /log, seed=,
;        printstep=, file=
;      Run a markov chain, optionally outputting to a file.  Both a functional
;      and procedural interface are provided.
;  mcmc::read_trials(file):  Read a set of trials from an output file.
;  mcmc::extract_stats(trialp, burning): Extract stats for each parameter from
;     a trials parameter array with the specified burn in.
;
; MODIFICATION HISTORY:
;  Created: 2005: Erin Sheldon, NYU.
;  Incorporated Hogg's more general method, and a bug fix thanks to Morad
;    Masjedi.  2006-10-14 Erin Sheldon, NYU
;
;-

FUNCTION mcmc::init
  return,1
END 

;docstart::mcmc::run
; NAME:
;  MCMC::RUN()
;  MCMC::RUN
;
; PURPOSE:
;  Run a Monte Carlo Markov Chain.  Both a functional and procedural interface
;  are provided.
;
; CALLING SEQUENCE:
;  trial_pars = mcmc->run(step_func, like_func, nstep, parguess, like=, 
;                        /log, seed=, file=)
;    OR
;  mcmc->run, step_func, like_func, nstep, parguess, trial_pars, like=, 
;        /log, seed=, file= 
;
; INPUTS:
;  step_func: A function that takes a step in the chain.
;  like_func: A function that computes the likelihood of an input set of
;             parameters. 
;  nstep: The number of steps to take in the Markov Chain.
;  parguess:  A starting guess at the parameters.
;
; KEYWORD PARAMETERS:
;  /log: "like_func" will calculate the log likelihood.  This is usually a
;        better way to go because of precision issues with exponentials in the
;        likelihood function.
;  seed: A starting seed.  If not input one is generated using IDL's builtin
;        seed generator from randomu.
;  printstep:  Print out a progress report every printstep trials.
;  file: Write the results to an output file rather than returning as a
;        variable. The value of trial_pars will be -1 in this case. The
;        results are read using the MCMC::READ_TRIALS() function.
;
; OUTPUTS:
;  trial_pars: An array [npars,ntrials] containing the trial steps through the
;              parameter space.  Will be -1 if file= is sent.
;
; OPTIONAL OUTPUTS:
;  like: The likelihood for each trial.
;
;
; EXAMPLES:
; 
; ; Fit a constant to n measurements.
;
; ; log likelihood.  
; function const_like, p
;   common mcmc_test_block, x, y, ivar, psigma, npars
;   chi2 = total( (y-p[0])^2*ivar )
;   return,-0.5*chi2
; end 
;
; function const_step, seed, pars
;   common mcmc_test_block, x, y, ivar, psigma, npars
;   return, pars + psigma*randomn(seed)
; end 
; function const_truepars, sigma=sigma
;   sigma = 1.0
;   return, 10.0
; end 
;
; ; For this test, just fake the data.
; pro const_setup
;   common mcmc_test_block, x, y, ivar, psigma, npars
;   ny = 10
;   pars = mcmc_test_const_truepars(sigma=sigma)
;   y = replicate(pars, ny)
;   yerr = replicate(sigma, ny)
;
;   ivar = 1.0/yerr^2
;   y[*] = y[*] + yerr*randomn(seed, ny)
;   psigma = yerr[0]
;   npars = 1
; end 
;
; IDL> mcmc=obj_new('mcmc')
; IDL> const_setup
; IDL> nstep = 10000
; IDL> parguess = 5.0 ; initial guess.  
; IDL> trials = mcmc->run('const_step', 'const_like', nstep, parguess, /log)
;
; ; Or you can write them to a file.  Better of it will take a long time
; ; or might crash.  This way you keep a record of what happened.
; IDL> mcmc->run, 'const_step', 'const_like', nstep, parguess, file=file, /log
; IDL> trials = mcmc->read_trials(file)
;
; MODIFICATION HISTORY:  
;  Created: 2005: Erin Sheldon, NYU.
;  Incorporated Hogg's more general method, and a bug fix thanks to Morad
;    Masjedi.  2006-10-14 Erin Sheldon, NYU
;
;docend::mcmc::run

function mcmc::run, step_func, like_func, nstep, parguess, like=like, log=log, seed=seed, file=file, printstep=printstep
  self->run, step_func, like_func, nstep, parguess, pars, like=like, log=log, seed=seed, file=file, printstep=printstep
  return, pars
end 
pro mcmc::run, step_func, like_func, nstep, parguess, pars, like=like, log=log, seed=seed, file=file, printstep=printstep

  if n_elements(seed) eq 0 then begin 
      seed = long( randomu(seed2)*10000 )
      tmp = randomu(seed)
  endif 

  npars   = n_elements(parguess)
  oldpars = float(parguess)

  ;; Should we write this to a file?
  if n_elements(file) ne 0 then begin 
      openw, lun, file, /get_lun

      ;; longs
      writeu, lun, long(nstep)
      writeu, lun, long(npars)

      pars = -1
  endif else begin 
      pars    = replicate(parguess[0],npars,nstep)
      like    = fltarr(nstep)
  endelse 

  for ii=0L,nstep-1L do begin

      self->_step, seed, oldpars, oldlike, step_func, like_func, newpars, newlike,$
        log=log

      if n_elements(file) ne 0 then begin 
          writeu, lun, float(newpars)
      endif else begin 
          pars[*,ii] = newpars
          like[ii] = newlike
      endelse 
      
      oldpars = newpars
      oldlike = newlike
      
      if n_elements(printstep) ne 0 then begin 
          if (ii mod  printstep) eq 0 then begin 
              print, $
                strjoin(strarr(21)+string(byte(8)),''), $
                'MCMC: ',100L*ii/nstep,' percent', $
                format= '($,A,A,I2.2,A)'
          endif 
      endif 
  endfor

  if n_elements(file) ne 0 then begin
      flush, lun
      free_lun, lun
  endif 

  if n_elements(printstep) ne 0 then begin 
      print, strjoin(strarr(21)+string(byte(8)),'')+'MCMC: done      '
  endif 
  return

end 

pro mcmc::_step, seed, pars, like, step_func, like_func, newpars, newlike, log=log

  if (not keyword_set(like)) then like= call_function(like_func,pars)

  newpars= call_function(step_func,seed,pars)
  newlike= call_function(like_func,newpars)

  if keyword_set(log) then likeratio= newlike-like $
  else likeratio= newlike/like

  randomnumber= randomu(seed)

  if keyword_set(log) then randomnumber= alog(randomnumber)
  if NOT ((newlike GT like) OR $
          (randomnumber LT likeratio)) then begin
      newpars= pars
      newlike= like
  endif

  return
end

;docstart::mcmc::read_trials
; NAME:
;  MCMC::READ_TRIALS()
;
; PURPOSE:
;  Read an output file from the MCMC ::run procedure.
;
; CALLING SEQUENCE:
;  trial_pars = mcmc->read_trials(file)
;
; INPUTS:
;  file: A file created by the MCMC::RUN procedure with file=file input.
;
; OUTPUTS:
;  trial_pars: An array [npars,ntrials] containing the trial steps through the
;              parameter space. 
;
;
;
; MODIFICATION HISTORY:  
;  Created: 2005: Erin Sheldon, NYU.
;
;docend::mcmc::read_trials


FUNCTION mcmc::read_trials, file
  openr, lun, file, /get_lun
  ntrial=0L
  npar=0L

  readu, lun, ntrial
  readu, lun, npar

  pars = fltarr(npar, ntrial)
  readu, lun, pars
  free_lun, lun

  return, pars

END 


;docstart::mcmc::extract_stats
; NAME:
;  MCMC::EXTRACT_STATS()
;
; PURPOSE:
;  Extract statistics from an MCMC chain with burn in.
;
; CALLING SEQUENCE:
;  stats = mcmc->extract_status(trial_pars, burnin)
;
; INPUTS:
;  trial_pars: The outputs of the MCMC chain.
;  burnin: Number of points from the beginning to skip for burn in.
;
; OUTPUTS:
;  Mean and standard deviation for each parameter.  A [npar,2] array.
;
; MODIFICATION HISTORY:  
;  Created: 2005: Erin Sheldon, NYU.
;
;docend::mcmc::extract_stats


FUNCTION mcmc::extract_stats, trialp, burnin

  IF n_elements(trialp) eq 0 or n_elements(burnin) EQ 0 THEN BEGIN 
      on_error, 2
      print,'-Syntax: parms = mcmc->extract_stats(trialp, burnin)'
      print
      message,'Halting'
  ENDIF 

  sst = size(trialp,/struct)

  if sst.n_dimensions ne 2 then begin 
      on_error,2
      message,'trialp must be a [npar,ntrial] array'
  endif 
  

  npar = sst.dimensions[0]
  ntrial = sst.dimensions[1]

  parms = dblarr(npar, 2)

  FOR i=0L, npar-1 DO BEGIN 
      mom = moment(trialp[i,burnin:ntrial-1], sdev=sd, /double)
      parms[i,0] = mom[0]
      parms[i,1] = sd
  ENDFOR 

  return, parms
END 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Some plotting routines
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Plot all the 2-d marginalizations

PRO mcmc::plot2d, burnin=burnin, mark=mark, nxbins=nxbins, nybins=nybins, names=names, _extra=_extra

  !p.multi=[0,2,2]


  IF n_elements(burnin) EQ 0 THEN burnin = 20
  IF n_elements(nxbins) EQ 0 THEN nxbins=30
  IF n_elements(nybins) EQ 0 THEN nybins=30

  IF n_elements(names) EQ 0 THEN BEGIN 

  ENDIF ELSE BEGIN
      
  ENDELSE 

  ntrial = self.ntrial

  parms = self->extract_stats(burnin)
  print,ntostr(parms[0,0])+' '+!plusminus+' '+ntostr(parms[0,1])
  print,ntostr(parms[1,0])+' '+!plusminus+' '+ntostr(parms[1,1])

  ;; 2d histogram
  histogram_2d, (*self.trial_struct).trialp[0,burnin:ntrial-1], (*self.trial_struct).trialp[1,burnin:ntrial-1], hist, $
    crap1, crap2, nxbins, nybins, /silent

  hx = total(hist.map, 2, /int)
  hy = total(hist.map, 1, /int)

  ;; Use our qmodel to get the model values
  crap = CALL_FUNCTION(self.qmodel, *self.x, *self.y, *self.ivar, reform(parms[*,0]), model=model)

  pplot, *self.x, *self.y, yerr=*self.yerr, psym=8, _extra=_extra
  pplot, *self.x, model, /overplot, color=!red

  plot, hist.xbins, hx, psym=10
  IF n_elements(mark) NE 0 THEN BEGIN 
      oplot, [ mark[0],mark[0] ], [0,1.e10], color=!red
  ENDIF 
  
  plot, hist.ybins, hy, psym=10
  IF n_elements(mark) NE 0 THEN BEGIN 
      oplot, [ mark[1],mark[1] ], [0,1.e10], color=!red
  ENDIF 


  Lrel = hist.map/max(hist.map)
  ;; convert chi^2 levels to likelihood ratio levels
  Llevels = exp( -0.5*!siglevels2 )

  image_contour, Lrel, hist.xbins, hist.ybins
  IF n_elements(mark) NE 0 THEN BEGIN 
      oplot, [ mark[0] ], [ mark[1] ], psym=7, symsize=2, color=!red
  ENDIF 

;  prob = hist.map/qgauss2d( hist.map, hist.xbins, hist.ybins, 100 )


  !p.multi=0

END 







;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; cleanup
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO mcmc__define
  struct = { $
             mcmc, $
             mcmc_dummy_var: 0 $
           }
END 
