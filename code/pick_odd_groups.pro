pro pick_odd_groups

dir="~alexie/Work/GroupCatalogs/"
group=mrdfits(dir+"cosmos_xgroups_20110209.fits",1)

minDist=50 ; kpc

d_MMGG200_as=distance(group.alpha_mmgg_scale,group.delta_mmgg_scale,group.alpha_mmgg_r200,group.delta_mmgg_r200)*3600.
d_MMGG200_kpc=d_MMGG200_as * (group.lensing_r200_mpc*1000./group.lensing_r200_as)

d_MLGGscale_as=distance(group.alpha_mmgg_scale,group.delta_mmgg_scale,group.alpha_mlgg_scale,group.delta_mlgg_scale)*3600.
d_MLGGscale_kpc=d_MLGGscale_as * (group.lensing_r200_mpc*1000./group.lensing_r200_as)

d_MLGG200_as=distance(group.alpha_mmgg_scale,group.delta_mmgg_scale,group.alpha_mlgg_r200,group.delta_mlgg_r200)*3600.
d_MLGG200_kpc=d_MLGG200_as * (group.lensing_r200_mpc*1000./group.lensing_r200_as)

sel_MMGG200=where(group.flag_include EQ 1 $
                  AND d_MMGG200_kpc GT minDist)

sel_MLGGscale=where(group.flag_include EQ 1 $
                  AND d_MLGGscale_kpc GT minDist)

sel_MLGG200=where(group.flag_include EQ 1 $
                  AND d_MLGG200_kpc GT minDist)

stop
end
