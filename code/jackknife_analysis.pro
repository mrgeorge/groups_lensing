pro analyze_file, file

readcol,file,run,nDiff,z,lx,mlx,mwl,sig_mwl,/silent
print,n_elements(run),n_elements(where(mwl LT mwl[0])),min(mwl),mean(mwl),stddev(mwl),mwl[0],nDiff[0]

end

pro jackknife_analysis,hiZ=hiZ

dir="../outfiles/bin_20_70_1000_7_emp_20110914/"

if(keyword_set(hiZ)) then zExt='_hiz' else zExt=''
bgg_scale_file=dir+"jackknife_results_20110914_mmgg_scale_bgg_scale"+zExt+".txt"
mmgg_r200_file=dir+"jackknife_results_20110914_mmgg_scale_mmgg_r200"+zExt+".txt"
bgg_r200_file=dir+"jackknife_results_20110914_mmgg_scale_bgg_r200"+zExt+".txt"

print,'mmgg_scale vs bgg_scale'
analyze_file,bgg_scale_file
print,'mmgg_scale vs mmgg_r200'
analyze_file,mmgg_r200_file
print,'mmgg_scale vs bgg_r200'
analyze_file,bgg_r200_file

stop
end
