function get_msun_lens, group, center

; return array of stellar masses for a given center type
; intended to be used by run_gg_offset to allow modeling of central
; point source component when lensing around different centers

mass=fltarr(n_elements(group))

case center of
   'mmgg_scale': mass=group.mmgg_scale_mstar
   'mmgg_r200': mass=group.mmgg_r200_mstar
   'mmgg2_r200': mass=group.mmgg2_r200_mstar
   'mlgg_scale': mass=group.mlgg_scale_mstar
   'bgg_scale': mass=group.mlgg_scale_mstar
   'mlgg_r200': mass=group.mlgg_r200_mstar
   'bgg_r200': mass=group.mlgg_r200_mstar
   'bcg1': mass=group.bcg_sm ; AL's old center
   'visual': mass=group.visual_mstar ; see mrg_center_notes_20110808_distilled.txt and weird_stack.pro
   'galaxy': mass=group.kevin_mstar
   else: mass[*]=0.              ; center is not a galaxy
endcase

return, mass
end
