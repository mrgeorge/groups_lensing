pro plot_chisq

fileDir="../outfiles/bin_20_70_1000_7_emp_20110914/"
;cenNames=['cn','cm','cf','xray','mmgg_scale','mmgg_r200','bgg_scale','bgg_r200']
cenNames=['mmgg_scale','bgg_scale']
lensOutFileArr=strcompress(fileDir+'center_'+cenNames+'.fits',/remove_all)

simpctable
plot,[30,1000],[0,1],/nodata,/xl,xst=1,color=!black,background=!white

for ii=0,n_elements(cenNames)-1 do begin
   str=mrdfits(lensOutFileArr[ii],1)
   get_ds_model,str.fit_type,str.p_mean,str.z_lens,str.msun_lens,str.plot_radius_kpc/1000.,cen_type='ps',/use_m200,star_term=star_term,nfw_term=nfw_term
   get_ds_model,str.fit_type2,str.p_mean2,str.z_lens,str.msun_lens,str.plot_radius_kpc/1000.,cen_type='ps',/use_m200,star_term=star_term,nfw_term=nfw_off

   if(str.msun_lens GT 0) then begin
      tot_cen=nfw_term+star_term
      tot_off=nfw_off+star_term
   endif else begin
      tot_cen=nfw_term
      tot_off=nfw_off
   endelse

   sel=where(str.plot_radius_kpc GT 30.)

   chi2_cen=(str.we1_mean-tot_cen)^2/str.we1_error^2
   chi2_off=(str.we1_mean-tot_off)^2/str.we1_error^2

   oplot,str.plot_radius_kpc[sel],chi2_cen[sel]/total(chi2_cen[sel]),color=ii
   oplot,str.plot_radius_kpc[sel],chi2_off[sel]/total(chi2_off[sel]),linestyle=2,color=ii
;   oplot,str.plot_radius_kpc[sel],tot_cen,color=ii
;   oplot,str.plot_radius_kpc[sel],tot_off,linestyle=2,color=ii
;   oploterror,str.plot_radius_kpc[sel],str.we1_mean,str.we1_error,ps=3,errcolor=ii

   print,cenNames[ii]
   print,chi2_cen[sel]
   stop
endfor

stop
end
