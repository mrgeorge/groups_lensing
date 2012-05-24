pro bootstrap_w, arr, weights, average, sigma, nSample
;  nSample=1000
  size=n_elements(arr)
  if(size LT nSample^(1./size)) then begin
     print, 'Too few elements in array for bootstrapping'
     stop
  endif

  sampleInd=round(randomu(seed,size,nSample)*(size-1))
;  avgArr=total(arr[sampleInd],1)/size
  avgArr=total(arr[sampleInd]*weights[sampleInd],1)/total(weights[sampleInd],1)
  average=mean(avgArr)
  sigma=stddev(avgArr)
end

pro centroid_uncertainties

; the code for getting centroid uncertainties was a bit buggy in
; find_bcg. This code will first check that the centroids agree with
; the find_bcg calculation, then recalculate the centroid
; uncertainties. Eventually it should be integrated back into
; find_bcg.pro for future catalogs.

dir="~/data/cosmos/code/"
groupFile=dir+"group5_20110914.fits"
groupOut=dir+"group5_cenunc_20110914_run20120518.fits"
acsFile=dir+"lensing18_20110914.fits"

group=mrdfits(groupFile,1)
group=group[where(group.flag_include EQ 1)]
acs=mrdfits(acsFile,1)

avg_alpha_cm=fltarr(n_elements(group))
avg_delta_cm=fltarr(n_elements(group))

for groupnum=0,n_elements(group)-1 do begin
   mem=where(acs.group_id_best EQ group[groupnum].id $
             AND acs.p_mem_best GT 0.5,nMem)

   if(nMem GT 0) then begin
      ; Unweighted member centroid - ALPHA_CN
      weight=replicate(1.,nMem)
      alpha_cn=total(acs[mem].alpha_j2000)/n_elements(mem)
      delta_cn=total(acs[mem].delta_j2000)/n_elements(mem)
;      alpha_cn_err=sqrt(total((acs[mem].alpha_j2000-alpha_cn)^2) / (nMem^2 - nMem))
;      delta_cn_err=sqrt(total((acs[mem].delta_j2000-delta_cn)^2) / (nMem^2 - nMem))
      bootstrap_w, acs[mem].alpha_j2000, weight, average, alpha_cn_err, 200
      bootstrap_w, acs[mem].delta_j2000, weight, average, delta_cn_err, 200
      group[groupnum].alpha_cn=alpha_cn
      group[groupnum].delta_cn=delta_cn
      group[groupnum].pos_err_cn=sqrt(alpha_cn_err*delta_cn_err)

      ; Stellar Mass weighted centroid - ALPHA_CM
      sm=where(acs[mem].kevin_mstar GT 2. AND acs[mem].kevin_mstar LT 15.,nsm)
      if(nsm GT 0) then begin
         weight=10.^(acs[mem[sm]].kevin_mstar)
         sumweight=total(weight)
         sumweight2=total(weight^2)
         alpha_cm=total(acs[mem[sm]].alpha_j2000 * weight)/total(weight)
         delta_cm=total(acs[mem[sm]].delta_j2000 * weight)/total(weight)
;         alpha_cm_err=sqrt(total(weight*(acs[mem[sm]].alpha_j2000-alpha_cm)^2)*sumweight2/(sumweight*(sumweight^2-sumweight2)))
;         delta_cm_err=sqrt(total(weight*(acs[mem[sm]].delta_j2000-delta_cm)^2)*sumweight2/(sumweight*(sumweight^2-sumweight2)))

;         alpha_cm_err=sqrt(total(weight^2*(acs[mem[sm]].alpha_j2000-alpha_cn)^2))/sumweight
;         delta_cm_err=sqrt(total(weight^2*(acs[mem[sm]].delta_j2000-delta_cn)^2))/sumweight
         
         bootstrap_w, acs[mem[sm]].alpha_j2000, weight, average_alpha_cm, alpha_cm_err, 200
         bootstrap_w, acs[mem[sm]].delta_j2000, weight, average_delta_cm, delta_cm_err, 200
         avg_alpha_cm[groupnum]=average_alpha_cm
         avg_delta_cm[groupnum]=average_delta_cm
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
      weight=10.^(-0.4*acs[mem].mag_auto)
      sumweight=total(weight)
      sumweight2=total(weight^2)
      alpha_cl=total(acs[mem].alpha_j2000 * weight)/total(weight)
      delta_cl=total(acs[mem].delta_j2000 * weight)/total(weight)
;      alpha_cl_err=sqrt(total(weight*(acs[mem].alpha_j2000-alpha_cl)^2)*sumweight2/(sumweight*(sumweight^2-sumweight2)))
;      delta_cl_err=sqrt(total(weight*(acs[mem].delta_j2000-delta_cl)^2)*sumweight2/(sumweight*(sumweight^2-sumweight2)))

      bootstrap_w, acs[mem].alpha_j2000, weight, average, alpha_cl_err, 200
      bootstrap_w, acs[mem].delta_j2000, weight, average, delta_cl_err, 200

      group[groupnum].alpha_cl=alpha_cl
      group[groupnum].delta_cl=delta_cl
      group[groupnum].pos_err_cl=sqrt(alpha_cl_err*delta_cl_err)
   endif
endfor

stop

mwrfits,group,groupOut,/create

sel=where(group.flag_include EQ 1)

degkpc=3600.*group[sel].lensing_r200_mpc*1000./group[sel].lensing_r200_as
print,'CN',median(group[sel].pos_err_cn*degkpc),mean(group[sel].pos_err_cn*degkpc)
print,'CM',median(group[sel].pos_err_cm*degkpc),mean(group[sel].pos_err_cm*degkpc)
print,'CL',median(group[sel].pos_err_cl*degkpc),mean(group[sel].pos_err_cl*degkpc)

stop
end
