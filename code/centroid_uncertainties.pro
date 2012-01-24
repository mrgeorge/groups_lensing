pro centroid_uncertainties

; the code for getting centroid uncertainties was a bit buggy in
; find_bcg. This code will first check that the centroids agree with
; the find_bcg calculation, then recalculate the centroid
; uncertainties. Eventually it should be integrated back into
; find_bcg.pro for future catalogs.

dir="~/data/cosmos/code/"
groupFile=dir+"group5_20110914.fits"
groupOut=dir+"group5_cenunc_20110914.fits"
acsFile=dir+"lensing18_20110914.fits"

group=mrdfits(groupFile,1)
acs=mrdfits(acsFile,1)

for groupnum=0,n_elements(group)-1 do begin
   mem=where(acs.group_id_best EQ group[groupnum].id $
             AND acs.p_mem_best GT 0.5,nMem)

   if(nMem GT 0) then begin
      ; Unweighted member centroid - ALPHA_CN
      alpha_cn=total(acs[mem].alpha_j2000)/n_elements(mem)
      delta_cn=total(acs[mem].delta_j2000)/n_elements(mem)
      alpha_cn_err=sqrt(total((acs[mem].alpha_j2000-alpha_cn)^2)/(n_elements(mem)-1))
      delta_cn_err=sqrt(total((acs[mem].delta_j2000-delta_cn)^2)/(n_elements(mem)-1))
      print,'CN difference:',distance(alpha_cn,delta_cn,group[groupnum].alpha_cn,group[groupnum].delta_cn)
      group[groupnum].alpha_cn=alpha_cn
      group[groupnum].delta_cn=delta_cn
      group[groupnum].pos_err_cn=sqrt(alpha_cn_err*delta_cn_err)

      ; Stellar Mass weighted centroid - ALPHA_CM
      sm=where(acs[mem].kevin_mstar GT 2. AND acs[mem].kevin_mstar LT 15.,nsm)
      if(nsm GT 0) then begin
         alpha_cm=total(acs[mem[sm]].alpha_j2000 * 10.^(acs[mem[sm]].kevin_mstar))/total(10.^(acs[mem[sm]].kevin_mstar))
         delta_cm=total(acs[mem[sm]].delta_j2000 * 10.^(acs[mem[sm]].kevin_mstar))/total(10.^(acs[mem[sm]].kevin_mstar))
         alpha_cm_err=sqrt(total(10^(acs[mem[sm]].kevin_mstar)*(acs[mem[sm]].alpha_j2000-alpha_cm)^2)*total(10^(acs[mem[sm]].kevin_mstar))/(total(10^(acs[mem[sm]].kevin_mstar))^2-total((10^(acs[mem[sm]].kevin_mstar))^2)))
         delta_cm_err=sqrt(total(10^(acs[mem[sm]].kevin_mstar)*(acs[mem[sm]].delta_j2000-delta_cm)^2)*total(10^(acs[mem[sm]].kevin_mstar))/(total(10^(acs[mem[sm]].kevin_mstar))^2-total((10^(acs[mem[sm]].kevin_mstar))^2)))
         print,'CM difference:',distance(alpha_cm,delta_cm,group[groupnum].alpha_cm,group[groupnum].delta_cm)
         group[groupnum].alpha_cm=alpha_cm
         group[groupnum].delta_cm=delta_cm
         group[groupnum].pos_err_cm=sqrt(alpha_cm_err*delta_cm_err)
      endif else begin
         print,'centroid_uncertainties: group['+string(groupnum,format='(I0)')+'] No objects with stellar mass in member selection.'
         group[groupnum].alpha_cm=-999
         group[groupnum].delta_cm=-999
         group[groupnum].pos_err_cm=-999
      endelse

      ; Luminosity weighted centroid - ALPHA_CL
      alpha_cl=total(acs[mem].alpha_j2000 * 10.^(-0.4*acs[mem].mag_auto))/total(10.^(-0.4*acs[mem].mag_auto))
      delta_cl=total(acs[mem].delta_j2000 * 10.^(-0.4*acs[mem].mag_auto))/total(10.^(-0.4*acs[mem].mag_auto))
      alpha_cl_err=sqrt(total(10^(-0.4*acs[mem].mag_auto)*(acs[mem].alpha_j2000-alpha_cl)^2)*total(10^(-0.4*acs[mem].mag_auto))/(total(10^(-0.4*acs[mem].mag_auto))^2-total((10^(-0.4*acs[mem].mag_auto))^2)))
      delta_cl_err=sqrt(total(10^(-0.4*acs[mem].mag_auto)*(acs[mem].delta_j2000-delta_cl)^2)*total(10^(-0.4*acs[mem].mag_auto))/(total(10^(-0.4*acs[mem].mag_auto))^2-total((10^(-0.4*acs[mem].mag_auto))^2)))
      print,'CL difference:',distance(alpha_cl,delta_cl,group[groupnum].alpha_cl,group[groupnum].delta_cl)
      group[groupnum].alpha_cl=alpha_cl
      group[groupnum].delta_cl=delta_cl
      group[groupnum].pos_err_cl=sqrt(alpha_cl_err*delta_cl_err)
   endif
endfor

mwrfits,group,groupOut,/create
end
