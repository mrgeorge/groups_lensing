pro single_model, strFile, chainFile, covPlotFile, $
                  noModel=noModel, $
                  ps=ps,sis=sis,tis=tis,rhotis=rhotis, $
                  off2dDelta=off2dDelta,off3dDelta=off3dDelta,off3dMax=off3dMax, $
                  mcfix=mcfix,mfix=mfix,cfix=cfix

if(keyword_set(mcfix)) then begin
   fitType = [$
             1,$                ; 0  M0    : baryonic mass
             2,$                ; 1  R_vir : NFW virial mass
             2,$                ; 2  C     : NFW concentration
             0,$                ; 3  alpha : fraction
             0,$                ; 4  bias
             0,$                ; 5  m_sigma
             1]                 ; 6  offset
endif else if(keyword_set(mfix)) then begin
   fitType = [$
             1,$                ; 0  M0    : baryonic mass
             2,$                ; 1  R_vir : NFW virial mass
             1,$                ; 2  C     : NFW concentration
             0,$                ; 3  alpha : fraction
             0,$                ; 4  bias
             0,$                ; 5  m_sigma
             1]                 ; 6  offset
endif else if(keyword_set(cfix)) then begin
   fitType = [$
             1,$                ; 0  M0    : baryonic mass
             1,$                ; 1  R_vir : NFW virial mass
             2,$                ; 2  C     : NFW concentration
             0,$                ; 3  alpha : fraction
             0,$                ; 4  bias
             0,$                ; 5  m_sigma
             1]                 ; 6  offset
endif else begin
   fitType = [$
             1,$                ; 0  M0    : baryonic mass
             1,$                ; 1  R_vir : NFW virial mass
             1,$                ; 2  C     : NFW concentration
             0,$                ; 3  alpha : fraction
             0,$                ; 4  bias
             0,$                ; 5  m_sigma
             1]                 ; 6  offset
endelse

if(n_elements(noModel) EQ 0) then begin
   run_ds_mcmc, strFile, fitType, rob_p_mean, rob_p_sigma, /slow, stackx=stackx, chainFile=chainFile,burnin=burnin,/noSave,ps=ps,sis=sis,tis=tis,rhotis=rhotis,off2dDelta=off2dDelta,off3dDelta=off3dDelta,off3dMax=off3dMax
endif

ds_cov_plots,chainFile,fitType,covPlotFile,burnin=burnin

end
