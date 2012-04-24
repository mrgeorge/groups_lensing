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

end
