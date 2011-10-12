function nfw_sigma_offset_function, x

common common_nfw_sigma_offset_function, p_f, zl_f, roff_f, mass_def_flag

if (mass_def_flag eq 0) then res=x*tabulate_nfw_sigma_offset(x, roff_f, p_f, zl_f,/r200)
if (mass_def_flag eq 1) then res=x*tabulate_nfw_sigma_offset(x, roff_f, p_f, zl_f,/r180)

return, res

end
