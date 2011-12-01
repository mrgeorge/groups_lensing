pro single_model, strFile, chainFile, covPlotFile, noModel=noModel

fitType = [$
1,$             ; 0  M0    : baryonic mass
1,$             ; 1  R_vir : NFW virial mass
1,$             ; 2  C     : NFW concentration
0,$             ; 3  alpha : fraction
0,$             ; 4  bias
0,$             ; 5  m_sigma
1]               ; 6  offset

if(n_elements(noModel) EQ 0) then $
   run_ds_mcmc, strFile, fitType, rob_p_mean, rob_p_sigma, /slow, stackx=stackx, chainFile=chainFile,burnin=burnin,/noSave $
else begin
   str=mrdfits(strFile,1)
   fitType=str.fit_type2
   rob_p_mean=str.p_mean2
   rob_p_sigma=str.p_sigma2
endelse

ds_cov_plots,chainFile,fitType,rob_p_mean,rob_p_sigma,covPlotFile,burnin=burnin

end
