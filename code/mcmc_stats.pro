pro mcmc_stats,pars,p_mean,p_sigma,p_mean_rob,p_sigma_rob,ml,ml_lo,ml_hi,nsig=nsig,niter=niter

; This is where robust mean and sigma are calculated (p_mean_rob, p_sigma_rob)
; The peak of the likelihood distribution (ml) is also calculated using a KDE

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

; Estimate maximum likelihood
; could just plot a histogram and find the peak, but let's try a KDE
; see https://en.wikipedia.org/wiki/Kernel_density_estimation
; uses a gaussian kernel with optimal bandwidth (assumes underlying
; distribution is also gaussian, but should generally work ok even if not)

ml=p_mean
ml_lo=ml
ml_hi=ml

nsample=sqrt(num)
for ii=0,np-1 do begin
   xarr=findgen(nsample)/(nsample-1)*(max(pars[ii,*])-min(pars[ii,*])) + min(pars[ii,*])
   yarr=replicate(0.,nsample)
   bandwidth=(4*p_sigma_rob[ii]^5/(3.*num))^(1./5)

   for jj=0L,num-1 do yarr+=exp(-(xarr-pars[ii,jj])^2/(2.*bandwidth^2))

   ml[ii]=xarr[where(yarr EQ max(yarr))]
   ml_lo[ii]=max(xarr[where(total(yarr,/cum)/total(yarr) LT 0.16)]) - ml[ii] ; 16th percentile (lower error bar for 68% range)
   ml_hi[ii]=min(xarr[where(total(yarr,/cum)/total(yarr) GT 0.84)]) - ml[ii] ; 84th percentile (upper error bar for 68% range)
endfor

return
end

