function get_ds_chisq,fit_type,p_mean,zl,msun,x,y,yerr,dof=dof,center=center,refcen=refcen,groupFile=groupFile,use_m200=use_m200,cen_type=cen_type,off_type=off_type
; Calculate chi^2
; currently just NFW + point source if point source is included in model
; dof keyword will be filled and returned if provided

; calculate expected model values at locations of data points
get_ds_model,fit_type,p_mean,zl,msun,x,cen_term=model_cen,nfw_term=model_nfw,$
             center=center,refcen=refcen,groupFile=groupFile,nfw_off=model_nfw_off,use_m200=use_m200,$
             cen_type=cen_type,off_type=off_type

if(keyword_set(center) AND keyword_set(refcen)) then $
   model_tot=model_nfw_off $
else model_tot=model_nfw

if(fit_type[0] NE 0) then model_tot+=model_cen

chisq=total((model_tot-y)^2/yerr^2)
dof=n_elements(x)-n_elements(where(fit_type EQ 1))

return, chisq
end
