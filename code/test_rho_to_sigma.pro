function interp_rho, r1, r2, rho2
; return rho2 sampled on the same log-spaced radii as r1
return, interpol(rho2,alog10(r2),alog10(r1))
end

function append_rho_nfw, rkpc, rho, rnfw_kpc, conc, zl, extFactor
; take the density profile from CONTRA and extend following NFW
; returns a density profile using the same spacing as rkpc (assumes
; log spacing), with a max radius = extFactor * max(rkpc)
; rkpc is modified in place

rbin=rkpc[1]/rkpc[0] ; assumes log spacing
oldMaxR=max(rkpc)
newMaxR=oldMaxR*extFactor
nNew=alog(newMaxR/rkpc[0])/alog(rbin)
nOld=n_elements(rkpc)

rNew=rkpc[0]*rbin^indgen(nNew)
rhoNew=nfw_rho(rNew/1000.,[rnfw_kpc/1000.,conc],zl) / 1.e9 ; NFW Msun/kpc^3

rNew[0:nOld-2]=rkpc[0:nOld-2] ; old radii aren't quite log-spaced so preserve them as is
rhoNew[0:nOld-2]=rho[0:nOld-2] ; keep old rho at the original radii

rkpc=rNew
return,rhoNew
end


pro run_contra, filename, MAC, conc, logSM, logMh, reff, rnfw

;MAC=0                         ; 0=Blumenthal, 1=Gnedin04
DM=1                          ; 1=NFW for initial halo profile
BAR=2                         ; 2=Hernquist (DeV) profile for stars
TRACE=0                       ; 2=Hernquist for vdisp tracer pop -> set to 0 to run faster
ANIS=0                        ; 0=isotropic
c=conc                        ; halo concentration
n_dm=0                        ; ignored based on DM!=2 and TRACE!=7
fb=10.^logSM/10.^logMh        ; baryon fraction within rvir
rb=reff/rnfw                  ; final baryon scalelength in units of rvir
n_b=0                         ; used for truncation in vdisp calculation
ra=0                          ; ignored based on ANIS=0

fmtI="(I3)"
fmtF="(F12.8)"

contraDir="contra/"

args=string(MAC,format=fmtI)+string(DM,format=fmtI)+string(BAR,format=fmtI)+string(TRACE,format=fmtI)+string(ANIS,format=fmtI)+string(conc,format=fmtF)+string(n_dm,format=fmtI)+string(fb,format=fmtF)+string(rb,format=fmtF)+string(n_b,format=fmtI)+string(ra,format=fmtI)

spawn,contraDir+"contra "+args+" > "+filename

end

pro get_halo_properties,logMh,zl,conc,rnfw_kpc
; halo size following get_ds_model
conc=get_zhao_conc(logMh,zl)
overdensity = get_overdensity( zl, /r200)
rho_crit = critical_density(zl)
factor   = alog10(double((4.0/3.0)*!Pi*overdensity*rho_crit))
r_log    = (1.0/3.0)*(logMh-factor)
rnfw_mpc     = 10.0^(r_log) ; Mpc
rnfw_kpc=rnfw_mpc*1000. ; kpc

end

pro get_rho,zl,conc,logSM,logMh,reff,rnfw_kpc,rkpc,rho1,rho2,rho3
; fill in rkpc, rho1, rho2, rho3 (Msun/kpc^3)

filename="contra.out"
MAC=0 ; Blumenthal
run_contra,filename,MAC,conc,logSM,logMh,reff,rnfw_kpc
readcol,filename,r1,r2,m1,m2,rho1,rho2,gam1,gam2,v_circ,sigma_los,comment='#'
MAC=1 ; Gnedin
run_contra,filename,MAC,conc,logSM,logMh,reff,rnfw_kpc
readcol,filename,r1,r3,m1,m3,rho1,rho3,gam1,gam3,v_circ,sigma_los,comment='#'

; convert densities to Msun / kpc^3
mnfw=10.^logMh ; Msun
rho1 *= mnfw/rnfw_kpc^3
rho2 *= mnfw/rnfw_kpc^3
rho3 *= mnfw/rnfw_kpc^3

; get densities sampled on the same radii
rho2=interp_rho(r1,r2,rho2) 
rho3=interp_rho(r1,r3,rho3) 

; extrapolate rho1 and rho2 beyond rvir following NFW
extFactor=10
rkpc=r1*rnfw_kpc ; kpc
rkpc_tmp=rkpc
rho1=append_rho_nfw(rkpc_tmp,rho1,rnfw_kpc,conc,zl,extFactor)
rkpc_tmp=rkpc
rho2=append_rho_nfw(rkpc_tmp,rho2,rnfw_kpc,conc,zl,extFactor)
rho3=append_rho_nfw(rkpc,rho3,rnfw_kpc,conc,zl,extFactor)

end

function hernquist, rkpc, reff, logSM
aa=reff/(1.+sqrt(2.)) ; note, this assumes reff is the 3d 1/2 mass radius
sm=10.^logSM ; Msun
rho=(sm/(2.*!pi)) * (aa/rkpc) / (rkpc + aa)^3 ; Msun / kpc^3
return, rho
end

pro get_profiles,zl,conc,logSM,logMh,reff,rnfw_kpc,rkpc,rho0,rho1,rho2,rho3,sigma0,sigma1,sigma2,sigma3,ds0,ds1,ds2,ds3
; fill in rkpc plus rho, sigma, ds profiles (0=gal, 1=uncontracted NFW, 2=Blumenthal, 3=Gnedin)

; produce and read in profiles from CONTRA
get_rho,zl,conc,logSM,logMh,reff,rnfw_kpc,rkpc,rho1,rho2,rho3

; galaxy
rho0=hernquist(rkpc,reff,logSM)
sigma0=rho_to_sigma(rho0,rkpc,rkpc)/1.e6 ; Msun / pc^2
ds0=sigma_to_ds(rkpc,sigma0)/1.e6 ; Msun / pc^2

; uncontracted
sigma1=rho_to_sigma(rho1, rkpc, rkpc)/1.e6 ; Msun / pc^2
ds1=sigma_to_ds(rkpc,sigma1)/1.e6 ; Msun / pc^2

; Blumenthal
sigma2=rho_to_sigma(rho2, rkpc, rkpc)/1.e6 ; Msun / pc^2 - note, using same radii as the non-contracted profile
ds2=sigma_to_ds(rkpc,sigma2)/1.e6 ; Msun / pc^2

; Gnedin
sigma3=rho_to_sigma(rho3, rkpc, rkpc)/1.e6 ; Msun / pc^2 - note, using same radii as the non-contracted profile
ds3=sigma_to_ds(rkpc,sigma3)/1.e6 ; Msun / pc^2

end

pro test_rho_to_sigma

; define some halo properties
zl=0.
logSM=11.
reff=10. ; kpc
logMh=12. ; log M200
get_halo_properties,logMh,zl,conc,rnfw_kpc ; fill in conc,rnfw_kpc

; get rho,sigma,ds for 0=gal, 1=uncontracted NFW, 2=Blumenthal, 3=Gnedin
get_profiles,zl,conc,logSM,logMh,reff,rnfw_kpc,rkpc,rho0,rho1,rho2,rho3,sigma0,sigma1,sigma2,sigma3,ds0,ds1,ds2,ds3


simpctable
!p.thick=3
!x.thick=3
!y.thick=3
!p.charsize=1.2
!p.charthick=1.2
!p.font=0
!p.multi=[0,3,2]

set_plot,'ps'
device,filename="lensing_ac_"+string(reff,format='(I0)')+"_"+string(logSM,format='(F4.1)')+"_"+string(logMh,format='(F4.1)')+".eps",/encapsul,/helvetica,/color, xsize=8,ysize=5,/inches

xr=[1,1000]

; plot all profiles
plot,rkpc,rho1,/xl,/yl,xr=xr,xtit='r (kpc)',ytit=textoidl("\rho [M_{sun} / kpc^3]")
oplot,rkpc,rho2,color=!blue
oplot,rkpc,rho3,color=!green
oplot,rkpc,rho0,color=!orange,linestyle=2

plot,rkpc,sigma1,/xl,/yl,xr=xr,xtit='R (kpc)',ytit=textoidl("\Sigma [M_{sun} / pc^2]")
oplot,rkpc,sigma2,color=!blue
oplot,rkpc,sigma3,color=!green
oplot,rkpc,sigma0,color=!orange,linestyle=2

plot,rkpc,ds1,/xl,/yl,xr=xr,xtit='R (kpc)',ytit=textoidl("\Delta\Sigma [M_{sun} / pc^2]")
oplot,rkpc,ds2,color=!blue
oplot,rkpc,ds3,color=!green
oplot,rkpc,ds0,color=!orange,linestyle=2
fit_type=[2,0,0,0,0,0,0]
get_ds_model,fit_type,[0],0,logSM,rkpc/1000,cen_type="ps",cen_term=cen_term
oplot,rkpc,cen_term,color=!magenta,linestyle=1

; plot diff between contracted profiles to uncontracted
plot,rkpc,(rho2-rho1)/rho1,/xl,xr=xr,xtit='r (kpc)',ytit=textoidl("(\rho_{AC} - \rho_{NFW}) / \rho_{NFW}")
oplot,rkpc,(rho2-rho1)/rho1,color=!blue
oplot,rkpc,(rho3-rho1)/rho1,color=!green
oplot,replicate(reff,2),!y.crange,color=!red,linestyle=1
oplot,replicate(rnfw_kpc,2),!y.crange,color=!red,linestyle=2

plot,rkpc,(sigma2-sigma1)/sigma1,/xl,xr=xr,xtit='R (kpc)',ytit=textoidl("(\Sigma_{AC} - \Sigma_{NFW}) / \Sigma_{NFW}")
oplot,rkpc,(sigma2-sigma1)/sigma1,color=!blue
oplot,rkpc,(sigma3-sigma1)/sigma1,color=!green
oplot,replicate(reff,2),!y.crange,color=!red,linestyle=1
oplot,replicate(rnfw_kpc,2),!y.crange,color=!red,linestyle=2

plot,rkpc,(ds2-ds1)/ds1,/xl,xr=[1,1000],xtit='R (kpc)',ytit=textoidl("(\Delta\Sigma_{AC} - \Delta\Sigma_{NFW}) / \Delta\Sigma_{NFW}")
oplot,rkpc,(ds2-ds1)/ds1,color=!blue
oplot,rkpc,(ds3-ds1)/ds1,color=!green
oplot,replicate(reff,2),!y.crange,color=!red,linestyle=1
oplot,replicate(rnfw_kpc,2),!y.crange,color=!red,linestyle=2

device,/close

stop

end

