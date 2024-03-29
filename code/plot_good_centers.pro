pro plot_good_centers, runStr

; this code is to plot the lensing signal calculated and used in model_good_centers
; while the MCMC chain is still running

plotDir="~/data/cosmos/groups_lensing/plots/bin_20_70_1000_7_emp_conc_cen/"
singlePlotFile=plotDir+"good_centers_lnl_rhotis_3dM"+runStr+".eps"

outDir="~/data/cosmos/groups_lensing/outfiles/bin_20_70_1000_7_emp_conc_cen/"
chainFile=outDir+"good_centers_lnl"+runStr+".chain"
lensOutFile=outDir+"good_centers_lnl"+runStr+".fits"

fitType = [$
          1,$                   ; 0  M0    : baryonic mass
          1,$                   ; 1  R_vir : NFW virial mass
          1,$                   ; 2  C     : NFW concentration
          0,$                   ; 3  alpha : fraction
          0,$                   ; 4  bias
          0,$                   ; 5  m_sigma
          1]                    ; 6  offset

; read chain file
mcmc=obj_new('mcmc')
pars=mcmc->read_trials(chainFile, like=lnl)
if(n_elements(burnin) EQ 0) then burnin=0.01 * n_elements(pars[0,*])
pars=pars[*,burnin:*]
lnl=lnl[burnin:*]

;; get rob_p_mean and rob_p_sigma from chain
;mcmc_stats,pars,p_mean,p_sigma,rob_p_mean,rob_p_sigma,mpost,mpost_lo,mpost_hi

; get point in chain with maximum likelihood
mlpars=pars[*,(where(lnl EQ max(lnl)))[0]]

; plot signal and model
plot_lensing_results,lensOutFile,singlePlotFile,mlpars,fitType,'rhotis','max3d',/use_m200,/models
stop
end
