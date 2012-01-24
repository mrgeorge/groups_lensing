pro mcmc_stats,pars,p_mean,p_sigma,p_mean_rob,p_sigma_rob,nsig=nsig,niter=niter

; This is where robust mean and sigma are calculated (p_mean_rob, p_sigma_rob)

if n_params() eq 0 then begin
    print,'-syntax mcmc_stats,pars,p_mean,p_sigma,p_mean_rob,p_sigma_rob,nsig=nsig,niter=niter'
    return
endif

if n_elements(nsig) eq 0 then nsig=4     ; 4 sigma clipping
if n_elements(niter) eq 0 then niter=5   ; number of iterations

; pars         = Markov chain : [n_params, nsteps]
sz=size(pars,/dim)
np=sz[0]                        ; number of parameters
num=float(sz[1])                ; number of steps
p_mean=total(pars,2)/num        ; mean   

pp=p_mean#replicate(1.0,num)
p_sigma=sqrt(total((pars-pp)^2,2)/num)  ; sigma

p_mean_rob=p_mean
p_sigma_rob=p_mean

for i=0L, np-1 do begin
    p=reform(pars[i,*])
    sigma_clip,p,mm,ss,nsig=nsig,niter=niter
    p_mean_rob[i]=mm
    p_sigma_rob[i]=ss
endfor

return
end

