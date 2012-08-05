pro radio_stack

; Stack list of candidate groups from Castignani/Chiaberge detected
; from FRII sources and galaxy overdensities.

fileDir="../radiogroups/"
idstr='123'
inputFile=fileDir+"candidates"+idstr+".txt"
groupFile=fileDir+"group_radio"+idstr+".fits"
lensOutFile=fileDir+"center_radio"+idstr+".fits"
plotFile=fileDir+"center_radio"+idstr+".eps"
infile_source='/Users/alexie/Work/Weak_lensing/GG_cat_2006/gglensing_source_v1.7.fits' ; Using the new catalog (photoz version 1.7)

; convert input text file to groupFile
readcol,inputFile,id,ra,dec,z,zlo,zhi,format='I,D,D,F,F,F',comment='#'

group={id:0L, $
       alpha_j2000:0.D, $
       delta_j2000:0.D, $
       zphot:0., $
       box:0., $
                                ; the following are just included for
                                ; code to work and not really used
       kevin_mstar:11.2, $
       hl_kpc:0., $
       e1_r:0., $
       e2_r:0. $
      }
group=replicate(group,n_elements(id))
group.id=id
group.alpha_j2000=ra
group.delta_j2000=dec
group.zphot=z

mwrfits,group,groupFile,/create
assign_boxes,groupFile,groupFile

innerRadiusKpc=20.
secondRadiusKpc=70.
maxRadiusKpc=1000.
nRadiusBins=7
minLensZ=0.
maxLensZ=3.0
minLensMass=12.
maxLensMass=15.
box_factor=20.
zscheme=2

run_gg_offset, infile_source, groupFile, lensOutFile, innerRadiusKpc, secondRadiusKpc, maxRadiusKpc, nRadiusBins, minLensZ, maxLensZ, minLensMass, maxLensMass, box_factor,zscheme,/usespecz,center=cenName,/emp ; no stackx

fitType = [$
2,$             ; 1  M0    : baryonic mass
1,$             ; 1  R_vir : NFW virial mass
1,$             ; 2  C     : NFW concentration
0,$             ; 3  alpha : fraction
0,$             ; 4  bias
0,$             ; 5  m_sigma
0 ]             ; 6  offset
run_ds_mcmc, lensOutFile, fitType, rob_p_mean, rob_p_sigma,/fast,/ps

plot_lensing_results,lensOutFile,plotFile,rob_p_mean,fitType,'ps',/use_m200,/models

end
