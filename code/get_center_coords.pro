pro get_center_coords, group,name,ra,dec,sel

case name of
   'xray': begin
      ra=group.alpha_ellipse
      dec=group.delta_ellipse
   end
   'cn': begin
      ra=group.alpha_cn
      dec=group.delta_cn
   end
   'cm': begin
      ra=group.alpha_cm
      dec=group.delta_cm
   end
   'cl': begin
      ra=group.alpha_cl
      dec=group.delta_cl
   end
   'cn_red': begin
      ra=group.alpha_cn_red
      dec=group.delta_cn_red
   end
   'cm_red': begin
      ra=group.alpha_cm_red
      dec=group.delta_cm_red
   end
   'cl_red': begin
      ra=group.alpha_cl_red
      dec=group.delta_cl_red
   end
   'cm_iter': begin
      ra=group.alpha_cm_iter
      dec=group.delta_cm_iter
   end
   'mmgg_scale': begin
      ra=group.alpha_mmgg_scale
      dec=group.delta_mmgg_scale
   end
   'mlgg_scale': begin
      ra=group.alpha_mlgg_scale
      dec=group.delta_mlgg_scale
   end
   'mmgg_r200': begin
      ra=group.alpha_mmgg_r200
      dec=group.delta_mmgg_r200
   end
   'mmgg2_r200': begin
      ra=group.alpha_mmgg2_r200
      dec=group.delta_mmgg2_r200
   end
   'mlgg_r200': begin
      ra=group.alpha_mlgg_r200
      dec=group.delta_mlgg_r200
   end
   'bcg1': begin ; AL's old center
      ra=group.alpha_bcg1
      dec=group.delta_bcg1
   end
   'visual': begin ; see mrg_center_notes_20110808_distilled.txt and weird_stack.pro
      ra=group.alpha_visual
      dec=group.delta_visual
   end
   'e1': begin
      ra=group.alpha_e1
      dec=group.delta_e1
   end
   'e2': begin
      ra=group.alpha_e2
      dec=group.delta_e2
   end
   'e3': begin
      ra=group.alpha_e3
      dec=group.delta_e3
   end
   'e4': begin
      ra=group.alpha_e4
      dec=group.delta_e4
   end
   'e5': begin
      ra=group.alpha_e5
      dec=group.delta_e5
   end
   'e6': begin
      ra=group.alpha_e6
      dec=group.delta_e6
   end
   'e7': begin
      ra=group.alpha_e7
      dec=group.delta_e7
   end
   'e8': begin
      ra=group.alpha_e8
      dec=group.delta_e8
   end
   'e9': begin
      ra=group.alpha_e9
      dec=group.delta_e9
   end
   'e10': begin
      ra=group.alpha_e10
      dec=group.delta_e10
   end
   'mock': begin
      ra=group.alpha_j2000
      dec=group.delta_j2000
   end
   'af': begin
      ra=group.alpha_j2000
      dec=group.delta_j2000
   end
endcase

sel=where(ra GT -360. and dec GT -90.)

end
