pro run_gg_offset, infile_source, infile_lens, outfile, type, minRadiusKpc, maxRadiusKpc, nRadiusBins, array, box_factor,zscheme,$
            bcg=bcg,$
            cleanbcg=cleanbcg,$
            notbcg=notbcg,$
            random_kpc=random_kpc,$
            random_ang=random_ang,$
            field=field,$
            group=group,$
            offset=offset,$
	    specz=specz,$
            pbzk=pbzk,$
            combine_z=combine_z,$
            shape1=shape1,$
            shape2=shape2,$
            xgroups=xgroups,$
            zgroups=zgroups,$
            highz=highz,$
            ml=ml,$
            sig2=sig2,$
            sig3=sig3,$
            modify=modify,$
            usespecz=usespecz,$
            center=center,$
            refcen=refcen,$
            agn=agn,$
            stackx=stackx, $               ; Stack the X axis according to r/R200
            emp_var=emp_var

;----------------------------------------------------------------------
; Various tests 
; zscheme 0 -> 68 % cut + sigma cut=0.1
; zscheme 1 -> 68 % cut + sigma cut=0.2
; zscheme 2 -> remove double peak objects + 68 % cut + sigma cut=0.1
; zscheme 3 -> remove double peak objects + 68 % cut + sigma cut=0.2
;----------------------------------------------------------------------

; -- THINGS TO DO STILL -----------------------------------------------
; Rederive the shape noise with new catalog and G factor
;----------------------------------------------------------------------

;----------------------------------------------------------------------
; NAME:
;       run_gglensing
; PURPOSE:
;       Call gglensing program
; EXPLANATION:
;
; CALLING SEQUENCE:
;       
; INPUTS:       
;       None
;
; OPTIONAL KEYWORD INPUT:
;
; OUTPUTS:      
;
; SIDE EFFECTS:  
;
; PROCEDURE:
;
; RESTRICTIONS:
;      
; MODIFICATION HISTORY:
;       Alexie Leauthaud	LAM	01.Nov.2006
;-----------------------------------------------------------------------

;common test_bin, tt_bin

sm_max = 13.0
sm_min = 8.0

; This is where the lenses lie:
z_min  = 0.2
z_max  = 1.2    ; for gg lensing need this: feb 09 -> changed to 1.0 here ...

print,'Run gglensing ...'

struct_source = mrdfits(infile_source, 1, /silent)
struct_lens   = mrdfits(infile_lens, 1, /silent) 


if (type eq -5) then begin
    sel = where(struct_lens.chandra_agn eq 1)
    struct_lens(sel).zphot = struct_lens(sel).chandra_zeta
endif

if (zscheme eq 2 or zscheme eq 3) then begin
    print,'Number density before:',n_elements(struct_source.alpha_j2000)/(1.64*60.0*60.0)
    sel=where(struct_source.zp2 lt 0)  ;select sources without a secondary peak
    print,'Number density after:',n_elements(struct_source(sel).alpha_j2000)/(1.64*60.0*60.0)
    struct_source = struct_source(sel)
endif

; For the EX specz AGN : change to specz
if (type eq -3 or type eq -4) then begin
    sel=where(struct_lens.eze_agn_specz eq 1)
    struct_lens(sel).zphot = struct_lens(sel).eze_specz
endif

; Replace with zspec for SL:
;sel=where(struct_lens.cecile_sl eq 1)
;if (sel[0] ne -1) then begin
;    print,'Replacing the SL zphot with specz'
;    struct_lens(sel).zphot=struct_lens(sel).cecile_sl_z
;endif


;------------------------------------------------------------------------
; AGN: CHANGE Redshifts for AGN here and make a cut on AGN quality
; key: <agn>
;------------------------------------------------------------------------

if ( keyword_set(agn) ) then begin
    sel=where(struct_lens.xmm_good eq 1 or struct_lens.chandra_good eq 1)
    struct_lens = struct_lens(sel)
    print,'>> Number of AGN :',n_elements(sel)

    ; Replace with XMM redshift when avaliable:
    sel=where(struct_lens.xmm_good eq 1)
    struct_lens(sel).zphot = struct_lens(sel).xmm_zeta

    ; Replace with Chandra redshift when avaliable:
    sel=where(struct_lens.chandra_good eq 1)
    struct_lens(sel).zphot = struct_lens(sel).chandra_zeta
endif

;------------------------------------------------------------------------
; Center study 
; key: <typecenter>
;------------------------------------------------------------------------

if ( keyword_set(xgroups) and keyword_set(center) ) then begin
   ; FOR THE CONC STUDY
   if (center eq 'conc') then begin
      sel=where(struct_lens.flag_include eq 1)
      struct_lens = struct_lens(sel)
      struct_lens.alpha_j2000 = struct_lens.alpha_mmgg_scale
      struct_lens.delta_j2000 = struct_lens.delta_mmgg_scale
   endif else begin
      if (NOT(keyword_set(refcen))) then begin
      ; FOR FULL STACKS OF INDIVIDUAL CENTERS
         print, 'SINGLE STACK CENTER: ',center
         get_center_coords,struct_lens,center,ra,dec,sel ;sel is just where this center is positive i.e. exists
         struct_lens.alpha_j2000=ra
         struct_lens.delta_j2000=dec
         struct_lens=struct_lens[sel]
         sel=where(struct_lens.flag_include EQ 1)
         struct_lens=struct_lens[sel]
         print,'NUM',n_elements(sel)
      endif else begin
      ; FOR COMPARING DIFFERENT CENTERS
         print, 'COMPARING DIFFERENT CENTERS: stacking on ',center,' w/ offset model relative to ',refcen
         get_center_coords,struct_lens,center,ra1,dec1
         get_center_coords,struct_lens,refcen,ra2,dec2
         offset_as=distance(ra1,dec1,ra2,dec2)*3600.
         offset_kpc=offset_as*struct_lens.lensing_r200_mpc*1000./struct_lens.lensing_r200_as
         min_offset_kpc=50.
         max_offset_kpc=200.
         min_mass=13.2
         max_mass=13.8
         struct_lens.alpha_j2000=ra1
         struct_lens.delta_j2000=dec1
         sel=where(struct_lens.flag_include EQ 1 $
                   AND ra1 GT 0 $
                   AND ra2 GT 0 $
                   AND offset_kpc GT min_offset_kpc $
;                   AND offset_kpc LT max_offset_kpc $
;                   AND struct_lens.lensing_m200 GT min_mass $
;                   AND struct_lens.lensing_m200 LT max_mass $
                  )
         struct_lens=struct_lens[sel]
         print,'NUM',n_elements(sel)
      endelse
   endelse
endif

    ;------------------------------------------------------------------------
    ; Change to BCG RA and DEC
    ; key: <bcg>, <cleanbcg>
    ;------------------------------------------------------------------------

    ;if ( keyword_set(xgroups) and keyword_set(bcg) ) then begin
    ;    print,'>>>  changing coordinates here ...'
    ;    struct_lens.alpha_j2000 = struct_lens.alpha_bcg1
    ;    struct_lens.delta_j2000 = struct_lens.delta_bcg1
    ;    sel=where(struct_lens.bcg_flag eq 2 or struct_lens.bcg_flag eq 3)   ; those that have been classified as a good BCG
    ;    print,'GROUP cut on BCG : ',n_elements(sel)
    ;    struct_lens=struct_lens(sel)
    ;endif

    ;if ( keyword_set(xgroups) and keyword_set(cleanbcg) ) then begin
    ;    print,'>>>  changing coordinates here ...'
    ;    struct_lens.alpha_j2000 = struct_lens.alpha_bcg1
    ;    struct_lens.delta_j2000 = struct_lens.delta_bcg1
    ;    sel=where(struct_lens.bcg_flag eq 2)   ; those that have been classified as a good BCG
    ;    print,'GROUP cut on BCG : ',n_elements(sel)
    ;    struct_lens=struct_lens(sel)
    ;endif

;------------------------------------------------------------------------
; Significance of X-ray groups
; key: <sig2>
;------------------------------------------------------------------------

if ( keyword_set(xgroups) and keyword_set(sig2) ) then begin
    sel=where(struct_lens.significance ge 2)
    print,'Significance GT 2 : ',n_elements(sel)
    struct_lens = struct_lens(sel)
endif

;------------------------------------------------------------------------
; ZCOSMOS groups
; key: <zgroups>
;------------------------------------------------------------------------

if ( keyword_set(zgroups) ) then begin
   sel=where(struct_lens.zphot ge 0.5 and struct_lens.zphot lt 0.8 and struct_lens.grp gt 0 and struct_lens.obs_rich gt 3 and struct_lens.completeness gt 0.35)
    print,' number of groups: ',n_elements(sel)
    struct_lens = struct_lens(sel)
endif

;------------------------------------------------------------------------
; High redshift structures
; key: <highz>
;------------------------------------------------------------------------

if ( keyword_set(highz) ) then begin
        print,'High z : no cuts'
        sel=where(struct_lens.id eq 28)
        struct_lens(sel).alpha_j2000 = 149.4507
        struct_lens(sel).delta_j2000 = 1.6707089

        sel=where(struct_lens.id eq 44)
        struct_lens(sel).alpha_j2000 = 149.80547
        struct_lens(sel).delta_j2000 = 1.8780739

        sel=where(struct_lens.id eq 97)
        struct_lens(sel).alpha_j2000 = 149.86639
        struct_lens(sel).delta_j2000 = 2.0628704

        sel=where(struct_lens.id eq 132)
        struct_lens(sel).alpha_j2000 = 149.51581
        struct_lens(sel).delta_j2000 = 2.1751037

        sel=where(struct_lens.id eq 263)
        struct_lens(sel).alpha_j2000 = 150.31661 
        struct_lens(sel).delta_j2000 = 2.8130993

        sel=where(struct_lens.id eq 276)
        struct_lens(sel).alpha_j2000 = 149.98908
        struct_lens(sel).delta_j2000 = 2.6887464
endif

;------------------------------------------------------------------------
; PBZK
; PBZK : Change zphot if Pbzk because Kevin's is truncated 
; key: <pbzk>
;------------------------------------------------------------------------

if keyword_set(pbzk) then begin
    struct_lens.zphot = struct_lens.pbzk_zphot
    struct_lens.kevin_mstar = alog10(struct_lens.oliv_mstar)
endif

;------------------------------------------------------------------------
; Some redshift cuts
;------------------------------------------------------------------------

; Just a zphot cut
sel=where(struct_lens.zphot ge array[2] and struct_lens.zphot le array[3])
print,''
print,'>> Just zphot cut : ',n_elements(sel),' objects'
print,'Redshift cut : ',array[2],array[3]
print,''

; Cut the catalog here
struct_lens=struct_lens(sel)

print,'Number of  Lenses :', n_elements(struct_lens)

;------------------------------------------------------------------------
; LENS SELECTION
;------------------------------------------------------------------------

select_lens = select_gg_lens(struct_lens, type, array)

if keyword_set(field) then begin
    select_lens = select_gg_lens(struct_lens, type, array,/field)
endif

if keyword_set(group) then begin
    select_lens = select_gg_lens(struct_lens, type, array,/group)
endif

;------------------------------------------------------------------------
; SELECT THE LENSES
;------------------------------------------------------------------------

struct_lens = struct_lens(select_lens)

if keyword_set(modify) then begin
    sel=where(struct_lens.id ne 178 and struct_lens.id ne 110)
    struct_lens = struct_lens(sel)
    ;110 bad 178 bad
    ;77 95 99 110 127 133 167 178 282 287 295 308
    print,struct_lens.im_num
    print,''
    print,struct_lens.id
    print,''
endif

n_lens = n_elements(struct_lens) & n_source = n_elements(struct_source)

print,""
print, "Number of sources:", n_source, "  Number of lens:", n_lens
print,""

;------------------------------------------------------------------------
; LENS Parameters
;------------------------------------------------------------------------

x_lens      =  struct_lens.alpha_j2000  ;Ra 
y_lens      =  struct_lens.delta_j2000  ;Dec
z_lens      =  struct_lens.zphot        ;Z
box_lens    =  struct_lens.box     

if (NOT(keyword_set(zgroups)) AND keyword_set(specz)) then begin
    zspec_lens  =  struct_lens.specz 
endif

;------------------------------------------------------------------------
; SOURCE Parameters
;------------------------------------------------------------------------

ident_source    = struct_source.ident
x_source        = struct_source.alpha_j2000  ;Ra 
y_source        = struct_source.delta_j2000  ;Dec
z_source        = struct_source.zphot        ;Z
z_source_low_68 = struct_source.photz_low_68
e1_source       = struct_source.gamma1       ;e1
e2_source       = struct_source.gamma2       ;e2 
var_e1_source   = struct_source.var_e1
var_e2_source   = struct_source.var_e2
box_source      = struct_source.box
sn_source       = struct_source.sn
mag_source      = struct_source.mag_auto
z_source_non_comb = struct_source.photoz_non_comb
z_source_non_comb_low68 = struct_source.photoz_non_comb_low_68
z_source_non_comb_high68 = struct_source.photoz_non_comb_high_68
source_good_specz = struct_source.good_specz

if keyword_set(random_ang) then begin
    ; Put all sources at z=1
    z_source[*]=1.0
endif

;------------------------------------------------------------------------
; Random Points
;------------------------------------------------------------------------

;if keyword_set(random_kpc) then begin

;    random_number = 100000
;    print,'------- RANDOM NUMBER ----------'
;    print,random_number

;    ra_max=max(struct_lens.alpha_j2000)  &  ra_min=min(struct_lens.alpha_j2000)
;    dec_max=max(struct_lens.delta_j2000) &  dec_min=min(struct_lens.delta_j2000)

;    if keyword_set(random_ang) then begin
    ; Put all sources at z=1
    ; Put all lenses at z=0.5
;    z_lens[*]=0.5
;    endif

;    rand_array,ra_max,ra_min,dec_max,dec_min,z_lens,box_factor,random_number,x_lens_rand,y_lens_rand,z_lens_rand,box_lens_rand

;    x_lens   = x_lens_rand
;    y_lens   = y_lens_rand
;    z_lens   = z_lens_rand
;    box_lens = box_lens_rand
;
;    print,''
;    print,'*** Number of Random Points : ***',n_elements(x_lens)
;    print,''

    ; Fill in other structures
;    mag_lens=fltarr(random_number)
;    msun_lens=fltarr(random_number)
   
;endif

;------------------------------------------------------------------------
; OPTIONS 
;------------------------------------------------------------------------

if keyword_set(bcg) then begin
	print,'> USE BCG'
endif
if keyword_set(RANDOM_KPC) then begin
	print,'> RANDOM POINTS'
endif

;------------------------------------------------------------------------
; Log bins
;------------------------------------------------------------------------

; Build radius array - radius_kpc[i] is the outer radius of the ith bin
;                      (zero is assumed for the inner radius of bin 0)
logFactor=(maxRadiusKpc/minRadiusKpc)^(1./(nRadiusBins-1))
radius_kpc=dblarr(nRadiusBins)
radius_kpc[0]=minRadiusKpc
for i=1,nRadiusBins-1 do begin
   radius_kpc[i]=radius_kpc[i-1]*logFactor
endfor

if keyword_set(stackx) then begin
    print,'> USING THE STACKX OPTION :'
    radius_kpc = radius_kpc/500.0
    print,radius_kpc
endif

;------------------------------------------------------------------------
; Call gglensing
;------------------------------------------------------------------------

; change coordinates here:
if keyword_set(highz) then begin  
   ;print,''
   ;print, struct_lens[0].alpha_j2000
   ;print, struct_lens[0].delta_j2000
   ;print, '150.23711  2.68332'
   ;print, struct_lens[1].alpha_j2000
   ;print, struct_lens[1].delta_j2000
   ;print,'150.28313  2.092607'

    ;if(type eq 999) then begin
    ;struct_lens[0].alpha_j2000 = 150.23711 
    ;struct_lens[0].delta_j2000 = 2.68332
    ;endif
    ;if(type eq 998) then begin
    ;struct_lens[0].alpha_j2000 = 150.28313
    ;struct_lens[0].delta_j2000 = 2.092607
    ;endif
    ;if(type eq 997) then begin
    ;struct_lens[0].alpha_j2000 = 150.23711 
    ;struct_lens[0].delta_j2000 = 2.68332

    ;struct_lens[1].alpha_j2000 = 150.28313
    ;struct_lens[1].delta_j2000 = 2.092607
    ;endif
endif

temp=create_struct($
'mag_lens',0D,$
'msun_lens',0D,$
'hl_kpc',0D,$
'x_lens',0D,$
'y_lens',0D,$
'z_lens',0D,$
'e1_lens',0D,$
'e2_lens',0D,$
'box_lens',0)

if keyword_set(xgroups) then begin
    temp=create_struct($
         'mag_lens',0D,$
         'msun_lens',0D,$
         'hl_kpc',0D,$
         'x_lens',0D,$
         'y_lens',0D,$
         'z_lens',0D,$
         'id',0,$
         'lx_scale',0.0,$
         'lensing_r200_mpc',0.0,$
         'e1_lens',0D,$
         'e2_lens',0D,$
         'box_lens',0)
endif

str_lens = replicate(temp,n_elements(x_lens))

temp=create_struct($
'ident_source',0L,$
'x_source',0D,$
'y_source',0D,$
'z_source',0D,$
'z_source_low_68',0D,$
'z_source_non_comb',0.0,$
'z_source_non_comb_low68',0.0,$
'z_source_non_comb_high68',0.0,$
'source_good_specz',0D,$
'sn_source',0.0,$
'mag_source',0.0,$
'box_source',0,$
'e1_source',0D,$
'e2_source',0D,$
'var_e1_source',0D,$
'var_e2_source',0D)

str_source = replicate(temp,n_elements(x_source))

if ((NOT keyword_set(zgroups)) AND (NOT keyword_set(highz)) ) then begin
    str_lens.x_lens           = x_lens
    str_lens.y_lens           = y_lens
    str_lens.z_lens           = z_lens
    str_lens.box_lens         = box_lens
    if NOT keyword_set(xgroups) then begin
        str_lens.hl_kpc       = struct_lens.hl_kpc  
        str_lens.msun_lens    = struct_lens.kevin_mstar
        str_lens.e1_lens      = struct_lens.e1_r    ;e1  (DONT USE GAMMA FOR THE LENSES)
        str_lens.e2_lens      = struct_lens.e2_r    ;e2  (REALLY JUST USING THE SHAPES HERE ...!)
    endif else begin
        str_lens.lx_scale     = struct_lens.lx_scale
        str_lens.lensing_r200_mpc = struct_lens.lensing_r200_mpc
        str_lens.id           = struct_lens.id
;        str_lens.msun_lens    = struct_lens.kevin_mstar
;        str_lens.msun_lens    = struct_lens.mmgg_scale_mstar
        str_lens.msun_lens    = get_msun_lens(struct_lens,center)
    endelse
    if keyword_set(specz) then begin
        str_lens.z_lens       = zspec_lens
    endif
endif else begin
    str_lens.x_lens           = x_lens
    str_lens.y_lens           = y_lens
    str_lens.z_lens           = z_lens
    str_lens.box_lens         = box_lens
endelse

str_source.ident_source             = ident_source
str_source.x_source                 = x_source
str_source.y_source                 = y_source
str_source.z_source                 = z_source
str_source.z_source_low_68          = z_source_low_68
str_source.z_source_non_comb        = z_source_non_comb
str_source.z_source_non_comb_low68  = z_source_non_comb_low68
str_source.z_source_non_comb_high68 = z_source_non_comb_high68
str_source.sn_source                = sn_source
str_source.mag_source               = mag_source
str_source.box_source               = box_source
str_source.e1_source                = e1_source
str_source.e2_source                = e2_source
str_source.var_e1_source            = var_e1_source
str_source.var_e2_source            = var_e2_source
str_source.source_good_specz        = source_good_specz

if keyword_set(offset) then begin
	print,'> OFFSET'
stop ; check that we're not confused with other offset code        
        delta=16     ; arcsec
        str_lens.x_lens =  str_lens.x_lens+delta
        str_lens.y_lens =  str_lens.y_lens+delta
endif

;--- Make some memory space ----
delvarx, struct_source, /FREE_MEM
delvarx, struct_lens,   /FREE_MEM

;----- CALL GG LENSING --------
if (keyword_set(shape1) or keyword_set(shape2)) then begin
    if keyword_set(shape1) then begin
        PRINT,'RUNNING GGLENSING_V6'
        gglensing_v6, str_lens, str_source, outfile, radius_kpc, box_factor, zscheme,/shape1
    endif
    if keyword_set(shape2) then begin
        PRINT,'RUNNING GGLENSING_V6'
        gglensing_v6, str_lens, str_source, outfile, radius_kpc, box_factor, zscheme,/shape2 
    endif
endif else begin
    if keyword_set(ml) then begin
        print,''
        PRINT,'RUNNING MAX LIKELIHOOD GGLENSING_V4'
        print,''
        ; This is r_min and r_max
        ; r_min should probably be set the minimum distance at which I think shape
        ; measurements are reliable
        if (bintype eq 200) then radius_kpc = [100,2000]
        if (bintype eq 201) then radius_kpc = [100,1000]
        if (bintype eq 202) then radius_kpc = [100,1500]
        if (bintype eq 203) then radius_kpc = [30,1500]
        gglensing_v4_ml, str_lens, str_source, outfile, radius_kpc, box_factor  
    endif else begin
        print,''
        PRINT,'RUNNING GGLENSING_V6'
       gglensing_v6, str_lens, str_source, outfile, radius_kpc, box_factor, zscheme,usespecz=usespecz,stackx=stackx,emp_var=keyword_set(emp_var)
    endelse
endelse
end
