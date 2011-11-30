pro run_gg_offset, infile_source, infile_lens, outfile, innerRadiusKpc, secondRadiusKpc, maxRadiusKpc, nRadiusBins, minLensZ, maxLensZ, minLensMass, maxLensMass, box_factor,zscheme,$
            xgroups=xgroups,$
            zgroups=zgroups,$
            usespecz=usespecz,$
            center=center,$
            refcen=refcen,$
            stackx=stackx, $               ; Stack the X axis according to r/R200
            emp_var=emp_var, $
            subhalo=subhalo                ; set msun_lens to be a fixed multiple of the stellar mass to test contribution of subhalo

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
;       removed lots of extraneous features
;        to streamline for the group centering study -MRG- 05 July 2011
;-----------------------------------------------------------------------

print,'Run gglensing ...'

struct_source = mrdfits(infile_source, 1, /silent)
struct_lens   = mrdfits(infile_lens, 1, /silent) 

if (zscheme eq 2 or zscheme eq 3) then begin
    print,'Number density before:',n_elements(struct_source.alpha_j2000)/(1.64*60.0*60.0)
    sel=where(struct_source.zp2 lt 0)  ;select sources without a secondary peak
    print,'Number density after:',n_elements(struct_source(sel).alpha_j2000)/(1.64*60.0*60.0)
    struct_source = struct_source(sel)
endif

;------------------------------------------------------------------------
; Center study 
;------------------------------------------------------------------------

if ( keyword_set(xgroups) and keyword_set(center) ) then begin
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
      struct_lens.alpha_j2000=ra1
      struct_lens.delta_j2000=dec1
      sel=where(struct_lens.flag_include EQ 1 $
                AND ra1 GT 0 $
                AND ra2 GT 0 $
                AND offset_kpc GT min_offset_kpc $
               )
      struct_lens=struct_lens[sel]
      print,'NUM',n_elements(sel)
   endelse
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
; LENS REDSHIFT SELECTION
;------------------------------------------------------------------------

; Just a zphot cut
sel=where(struct_lens.zphot GE minLensZ $
          AND struct_lens.zphot LT maxLensZ)
print,''
print,'>> Just zphot cut : ',n_elements(sel),' objects'
print,'Redshift cut : ',minLensZ,maxLensZ
print,''

; Cut the catalog here
struct_lens=struct_lens(sel)

print,'Number of  Lenses :', n_elements(struct_lens)

;------------------------------------------------------------------------
; LENS MASS SELECTION
;------------------------------------------------------------------------

zphot=struct_lens.zphot
lx=struct_lens.lx_scale

select_lens=where(struct_lens.lensing_m200 GE minLensMass $
                  AND struct_lens.lensing_m200 LE maxLensMass)
print,'After lens mass cut: ',minLensMass,maxLensMass

struct_lens = struct_lens(select_lens)

print,'Number of lenses :',n_elements(struct_lens)

n_lens=n_elements(struct_lens)
n_source=n_elements(struct_source)

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

;------------------------------------------------------------------------
; Log bins
;------------------------------------------------------------------------

; Build radius array - radius_kpc[i] is the outer radius of the ith bin
;                      (zero is assumed for the inner radius of bin 0)
logFactor=(maxRadiusKpc/secondRadiusKpc)^(1./(nRadiusBins-2))
radius_kpc=dblarr(nRadiusBins)
radius_kpc[0]=innerRadiusKpc
radius_kpc[1]=secondRadiusKpc
radius_kpc[2:nRadiusBins-1]=radius_kpc[1]*logFactor^(indgen(nRadiusBins-2)+1)
for i=2,nRadiusBins-1 do begin
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
   if(n_elements(subhalo) EQ 0) then $
      str_lens.msun_lens    = get_msun_lens(struct_lens,center) $
   else $
      str_lens.msun_lens    = alog10(subhalo) + get_msun_lens(struct_lens,center) ; for miscentering tests, in case candidate central is actually a satellite with a massive subhalo, to see how it affects the halo fit.
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

;--- Make some memory space ----
delvarx, struct_source, /FREE_MEM
delvarx, struct_lens,   /FREE_MEM

;----- CALL GG LENSING --------
print,''
PRINT,'RUNNING GGLENSING_V6'
gglensing_v6, str_lens, str_source, outfile, radius_kpc, box_factor, zscheme,usespecz=usespecz,stackx=stackx,emp_var=keyword_set(emp_var)

end
