pro runS82

dir="s82/"
nFiles=4
strFiles=dir+'test'+string(indgen(nFiles)+100,format='(I03)')+'.fits'
chainFiles=dir+'test'+string(indgen(nFiles)+100,format='(I03)')+'.chain'
plotFiles=dir+'test'+string(indgen(nFiles)+100,format='(I03)')+'.eps'

fitType = [$
          1,$                   ; 0  M0    : baryonic mass
          1,$                   ; 1  R_vir : NFW virial mass
          1,$                   ; 2  C     : NFW concentration
          0,$                   ; 3  alpha : fraction
          0,$                   ; 4  bias
          0,$                   ; 5  m_sigma
          1]                    ; 6  offset

for ii=0,n_elements(strFiles)-1 do begin
   run_ds_mcmc, strFiles[ii], fitType, rob_p_mean, rob_p_sigma, /fast, chainFile=chainFiles[ii],burnin=burnin,/rhotis,/off3dMax

   plot_lensing_results, strFiles[ii], plotFiles[ii], rob_p_mean, fitType, /use_m200, /models
endfor

end
