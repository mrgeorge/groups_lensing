pro model_good_centers

; latest group file
dir="~/data/cosmos/code/"
groupFile=dir+"group4_20110914.fits"
group=mrdfits(groupFile,1)

; copy box numbers from previous version
oldFile='/Users/alexie/Work/GroupCatalogs/cosmos_xgroups_20110209.fits' ; group catalog with centers
old=mrdfits(oldFile,1)

group.box=old.box

; select groups with galaxy centers that agree
sel=where(group.id_mmgg_scale EQ group.id_mlgg_scale $
          AND group.id_mmgg_scale EQ group.id_mmgg_r200 $
          AND group.id_mmgg_scale EQ group.id_mlgg_r200 $
          AND group.flag_include EQ 1, nSel)

print,'N groups',nSel

goodFile=dir+"good_centers_lnl_20110914.fits"
mwrfits,group[sel],goodFile,/create


; measure lensing signal for these groups

infile_source='/Users/alexie/Work/Weak_lensing/GG_cat_2006/gglensing_source_v1.7.fits' ; Using the new catalog (photoz version 1.7)
outDir="~/data/cosmos/groups_lensing/outfiles/bin_20_70_1000_7_emp_conc_cen/"
lensOutFile=outDir+"good_centers_lnl.fits"
innerRadiusKpc=20.
secondRadiusKpc=70.
maxRadiusKpc=1000.
nRadiusBins=7
minLensZ=0.
maxLensZ=1.0
minLensMass=12.
maxLensMass=15.
box_factor=20.
zscheme=2
cenName="mmgg_scale"

run_gg_offset, infile_source, goodFile, lensOutFile, innerRadiusKpc, secondRadiusKpc, maxRadiusKpc, nRadiusBins, minLensZ, maxLensZ, minLensMass, maxLensMass, box_factor,zscheme,/xgroups,/usespecz,center=cenName,/emp_var


; model the lensing signal

fitType = [$
          1,$                   ; 0  M0    : baryonic mass
          1,$                   ; 1  R_vir : NFW virial mass
          1,$                   ; 2  C     : NFW concentration
          0,$                   ; 3  alpha : fraction
          0,$                   ; 4  bias
          0,$                   ; 5  m_sigma
          1]                    ; 6  offset
chainFile=outDir+"good_centers_lnl.chain"
plotDir="~/data/cosmos/groups_lensing/plots/bin_20_70_1000_7_emp_conc_cen/"
covPlotFile=plotDir+"good_centers_lnl_rhotis_3dM.cov.eps"

run_ds_mcmc, lensOutFile, fitType, rob_p_mean, rob_p_sigma, /slow, chainFile=chainFile,burnin=burnin,/rhotis,/off3dMax

; make covar plot
ds_cov_plots,chainFile,fitType,covPlotFile,burnin=burnin

; plot lensing signal with model
;singlePlotFile=plotDir+"good_centers_rhotis_3dM.eps"
;plot_lensing_results,lensOutFile,singlePlotFile,rob_p_mean,fitType,/use_m200,/models,cen_type='rhotis',off_type='max3d'
plot_good_centers

end
