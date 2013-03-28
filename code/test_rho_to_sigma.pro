function interp_rho_f, r_i, r_f, rho_f
; return rho_f sampled on the same log-spaced radii as r_i
return, interpol(rho_f,alog10(r_f),alog10(r_i))
end

function append_rho_nfw, rkpc, rho, rnfw_mpc, conc, zl, extFactor
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
rhoNew=nfw_rho(rNew/1000.,[rnfw_mpc,conc],zl) / 1.e9 ; NFW Msun/kpc^3

rNew[0:nOld-2]=rkpc[0:nOld-2] ; old radii aren't quite log-spaced so preserve them as is
rhoNew[0:nOld-2]=rho[0:nOld-2] ; keep old rho at the original radii

rkpc=rNew
return,rhoNew
end


pro run_contra, filename, conc, logSM, logMh, reff, rnfw

MAC=0                         ; 0=Blumenthal, 1=Gnedin04
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

contraDir="~/sshfs/riemann/home/mgeorge/flow/code/contra/"

args=string(MAC,format=fmtI)+string(DM,format=fmtI)+string(BAR,format=fmtI)+string(TRACE,format=fmtI)+string(ANIS,format=fmtI)+string(conc,format=fmtF)+string(n_dm,format=fmtI)+string(fb,format=fmtF)+string(rb,format=fmtF)+string(n_b,format=fmtI)+string(ra,format=fmtI)

spawn,contraDir+"contra "+args+" > "+filename

end

pro test_rho_to_sigma

; define some halo properties
zl=0.
conc=3.66868
logMh=13. ; log M200
logSM=11.
reff=5. ; kpc

mnfw=10.^logMh ; Msun
; halo size following get_ds_model
overdensity = get_overdensity( zl, /r200)
rho_crit = critical_density(zl)
factor   = alog10(double((4.0/3.0)*!Pi*overdensity*rho_crit))
r_log    = (1.0/3.0)*(logMh-factor)
rnfw_mpc     = 10.0^(r_log) ; Mpc
rnfw_kpc=rnfw_mpc*1000. ; kpc


; produce and read in profiles from CONTRA
filename="contra.out"
run_contra,filename,conc,logSM,logMh,reff,rnfw_kpc
readcol,filename,r_i,r_f,m_i,m_f,rho_i,rho_f,gam_i,gam_f,v_circ,sigma_los,comment='#'


; convert densities to Msun / kpc^3
rho_i *= mnfw/rnfw_kpc^3
rho_f *= mnfw/rnfw_kpc^3

; get densities sampled on the same radii
rho_f=interp_rho_f(r_i,r_f,rho_f) 

; extrapolate rho_i and rho_f beyond rvir following NFW
rkpc=r_i*rnfw_kpc ; kpc
rkpc_tmp=rkpc
extFactor=10
rho_i=append_rho_nfw(rkpc_tmp,rho_i,rnfw_mpc,conc,zl,extFactor)
rho_f=append_rho_nfw(rkpc,rho_f,rnfw_mpc,conc,zl,extFactor)


; uncontracted
sigma1=rho_to_sigma(rho_i, rkpc, rkpc)/1.e6 ; Msun / pc^2

; contracted
sigma2=rho_to_sigma(rho_f, rkpc, rkpc)/1.e6 ; Msun / pc^2 - note, using same radii as the non-contracted profile

; NFW (uncontracted, calculated a different way to check for consistency)
sigma3=nfw_sigma(rkpc/1000.,[rnfw_mpc,conc],zl,/r200)


simpctable
!p.charsize=1.5

; compare contracted and uncontracted 3d density profiles
plot,rkpc,(rho_f-rho_i)/rho_i,/xl,ps=-4,xtit='r (kpc)',ytit=textoidl("(\rho_{AC} - \rho_{NFW}) / \rho_{NFW}")

; NFW (uncontracted, calculated a different way)
rho3=nfw_rho(rkpc/1000.,[rnfw_mpc,conc],zl) / 1.e9 ; Msun/kpc^3
oplot,rkpc,(rho_i - rho3)/rho3,color=!green
oplot,replicate(rnfw_kpc,2),!y.crange,linestyle=1

stop

; compare contracted and uncontracted 2d density profiles
plot,rkpc,(sigma2-sigma1)/sigma1,/xl,ps=-4,xtit='R (kpc)',ytit=textoidl("(\Sigma_{AC} - \Sigma_{NFW}) / \Sigma_{NFW}")
oplot,rkpc,(sigma3-sigma1)/sigma1,ps=-4,color=!green
oplot,replicate(rnfw_kpc,2),!y.crange,linestyle=1

stop






stop
end

