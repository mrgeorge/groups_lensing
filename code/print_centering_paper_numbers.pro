pro print_centering_paper_numbers
  
; print the number and density of source galaxies used in analysis for
; centering paper.
; code follows run_gg_offset.pro, which is called by full_stack and
; diff_stack

infile_source='/Users/alexie/Work/Weak_lensing/GG_cat_2006/gglensing_source_v1.7.fits' ; Using the new catalog (photoz version 1.7)

; zscheme 2 -> remove double peak objects + 68 % cut + sigma cut=0.1
zscheme=2

struct_source = mrdfits(infile_source, 1, /silent)

print,'N sources before cut',n_elements(struct_source)

if (zscheme eq 2 or zscheme eq 3) then begin
    print,'Number density before (#/sq arcmin):',n_elements(struct_source.alpha_j2000)/(1.64*60.0*60.0)
    sel=where(struct_source.zp2 lt 0)  ;select sources without a secondary peak
    print,'N sources after zp2 cut',n_elements(sel)
    print,'Number density after zp2 cut (#/sq arcmin):',n_elements(struct_source(sel).alpha_j2000)/(1.64*60.0*60.0)
    struct_source = struct_source(sel)
endif


; print number of group members, w/ and w/o speczs

dir="~/data/cosmos/code/"
group=mrdfits(dir+"group5_20110914.fits",1)
acs=mrdfits(dir+"lensing18_20110914.fits",1)
group=group[where(group.flag_include EQ 1)]

nphotmemtot=0
nphotmemspec=0
nspecmemtot=0
nspecmemspec=0
for ii=0,n_elements(group)-1 do begin
   ; photoz mem
   tmp=where(acs.p_mem_best GT 0.5 AND acs.group_id_best EQ group[ii].id,nphotmem)
   nphotmemtot+=nphotmem
   if(nphotmem GT 0) then begin
      tmp=where(acs[tmp].good_specz EQ 1,nmem)
      nphotmemspec+=nmem
   endif

   ; specz mem
   tmp=where(acs.p_mem_best_specz GT 0.5 AND acs.group_id_best_specz EQ group[ii].id,nspecmem)
   nspecmemtot+=nspecmem
   if(nspecmem GT 0) then begin
      tmp=where(acs[tmp].good_specz EQ 1,nmem)
      nspecmemspec+=nmem
   endif
endfor
print,'mean number of (photoz) members (group.n_mem)',mean(group.n_mem)
print,'mean number of (photoz) members (counting)',float(nphotmemtot)/n_elements(group)
print,'mean number of (photoz) members with specz (counting)',float(nphotmemspec)/n_elements(group)
print,'mean number of (specz) members (counting)',float(nspecmemtot)/n_elements(group)
print,'mean number of (specz) members with specz (counting)',float(nspecmemspec)/n_elements(group)


; mean MMGG_scale - X-ray offset
dist_mmgg_xray=distance(group.alpha_mmgg_scale,group.delta_mmgg_scale,group.alpha_ellipse,group.delta_ellipse)*group.lensing_r200_mpc/group.lensing_r200_as * (1000.*3600.)
print,'mean MMGG_scale - Xray offset (kpc)',mean(dist_mmgg_xray)
print,'median MMGG_scale - Xray offset (kpc)',median(dist_mmgg_xray)


; ML pars in cov fit
runStr="_07"
outDir="~/data/cosmos/groups_lensing/outfiles/bin_20_70_1000_7_emp_conc_cen/"
chainFile=outDir+"good_centers_lnl"+runStr+".chain"
mcmc=obj_new('mcmc')
pars=mcmc->read_trials(chainFile, like=lnl)
if(n_elements(burnin) EQ 0) then burnin=0.01 * n_elements(pars[0,*])
pars=pars[*,burnin:*]
lnl=lnl[burnin:*]
; get point in chain with maximum likelihood
mlpars=pars[*,(where(lnl EQ max(lnl)))[0]]
print,'ML pars (Msub, Mh, c, sigoff)',mlpars
zl=0.516
; following is taken from get_ds_model
M0=mlpars[0]
Mnfw=mlpars[1]
conc=mlpars[2]
offset=mlpars[3]
overdensity = get_overdensity( zl, /r200)
rho_crit = critical_density(zl)
factor   = alog10(double((4.0/3.0)*!Pi*overdensity*rho_crit))
r_log    = (1.0/3.0)*(Mnfw-factor)
rnfw     = 10.0^(r_log)
r_eff=0.005                  ; 5 kpc in Mpc
r_core=0.0001                   ; 0.1 kpc - to avoid singularity
minx=0.001
x_mpc_rho=minx*10.^(findgen(1000)/999.*3.3)
xsel=where(x_mpc_rho LT offset, nSel)
if(nSel GT 1) then begin
   ; define 3d densities for halo and subhalo
   halo_nfw_rho_off=nfw_rho(offset-x_mpc_rho[xsel],[rnfw,conc],zl) ; halo density along line connecting subhalo center with halo center, starting from the subhalo center

   rho0_pis=10.^(M0) / (4.*!pi*r_core^2 * (r_eff-r_core*atan(r_eff/r_core)))
   sub_pis_rho=rho0_pis / (1.+x_mpc_rho[xsel]^2/r_core^2)

   ; find truncation radius
   trunc_ind=min(where(halo_nfw_rho_off GT sub_pis_rho,nTrunc)) ; find where the halo starts to dominate the density
   if(nTrunc EQ 0) then begin                                   ; the subhalo dominates, don't truncate
      print,'GET_DS_MODEL: subhalo dominates halo, setting r_cut=rnfw'
      r_cut=rnfw                ; set r_cut = offset instead ?
   endif else if(trunc_ind EQ 0) then begin
                                ; halo dominates even at minx
      r_cut=minx
   endif else begin
                                ; interpolate over x_mpc_rho to find where densities are equal
      r_cut=(10.^(interpol(alog10(x_mpc_rho[xsel]),alog10(sub_pis_rho)-alog10(halo_nfw_rho_off),0.)))[0]
   endelse
endif
print, 'Truncation radius',r_cut
print, 'R200',rnfw
print, 'Rs',rnfw/conc
 

; mean centroid uncertainties
groupunc=mrdfits("~/data/cosmos/code/group5_cenunc_20110914_run20120518.fits",1)
degkpc=groupunc.lensing_r200_mpc*1000.*3600./group.lensing_r200_as
print,'Mean CN unc',mean(degkpc*groupunc.pos_err_cn)
print,'Mean CM unc',mean(degkpc*groupunc.pos_err_cm)
print,'Mean CL unc',mean(degkpc*groupunc.pos_err_cl)

end
