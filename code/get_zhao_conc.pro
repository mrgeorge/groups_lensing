function get_zhao_conc, mass, z

;>   m in log10 units
;>   For the mass 200*rho_crit
;>   Made by the program  make_conc_sav_file

common zhao_params, mc_mass_array, mc_conc_array

if NOT (keyword_set(mc_mass_array)) then begin
;    file='/Users/alexie/idl/Jeremy/Zhao/mc_mass_array.sav'
    file='zhao_mc_mass_array.sav'
    restore,filename=file
;    file='/Users/alexie/idl/Jeremy/Zhao/mc_conc_array.sav'
    file='zhao_mc_conc_array.sav'
    restore,filename=file
endif

zmin    =  0.0
zmax    =  1.5
zbin    =  0.02
nbin    =  floor(1+(zmax-zmin)/zbin)
zarray  =  findgen(nbin)*zbin
nn_z    =  n_elements(zarray)

; find the closest z: -> zindex
aa      = abs(zarray - z)
m       = min(aa,zindex)

mass_array = reform(mc_mass_array[zindex,*])

nn = n_elements(mass_array)

if ((mass le mass_array[0]) OR (mass ge mass_array[nn-1])) then begin
    print,'Mass out of range in get_zhao_conc.pro'
    print,'min = ',mass_array[0]
    print,'max = ',mass_array[nn-1]
    stop
endif

aa = abs(mass_array - mass)
m  = min(aa,mindex)

closest_mass = mass_array[mindex]
;print,'closest mass',closest_mass

if (abs(closest_mass - mass) lt 0.001) then begin
    res = mc_conc_array[zindex,mindex]
endif else begin
    if ((closest_mass - mass) gt 0) then index2 = mindex-1
    if ((closest_mass - mass) lt 0) then index2 = mindex+1

    ; Linear interpolation for Conc
    y1 = mc_conc_array(zindex, mindex)  ;! log value
    y2 = mc_conc_array(zindex, index2)
    x1 = mass_array[mindex]
    x2 = mass_array[index2]

    coef_a = (y1-y2)/(x1-x2)
    coef_b = y1-(coef_a*x1)

    res = (coef_a*mass)+coef_b
endelse

return,res

end
