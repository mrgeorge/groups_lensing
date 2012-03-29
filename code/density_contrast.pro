function density_contrast, z, bryan=bryan

;+
; NAME:
;   DENSITY_CONTRAST
;
; AUTHOR:
;   Alexie Leauthaud 2007
;
; PURPOSE:
;    Calculates density contrast, Delta_c, for given cosmology
;    Assumes flat universe. 
;
; INPUTS:
;    z     Redshift 
;		   
; KEYWORDS:
;    Bryan -> use Bryan and Norman 1998
;
; OUTPUT:
;    Virial Overdensity
;
; MODIFICATION HISTORY:
;
;-     

;-------------------------------------------------------------------------
; Based on <Pierpaoli, Scott & White 2000>
; Valid for 0.2 < Omega_m < 1.1 and 0 < Omega_l < 1
;-------------------------------------------------------------------------

; Get Omega_m at this redshift
Omega_m = get_omegam(z)

; Omega_l
Omega_l = 1-Omega_m

;-------------------------------------------------------------------------
; Pierpaoli
;-------------------------------------------------------------------------

x = Omega_m-0.2
y = Omega_l

A = fltarr(5,5)   ; i,j

A[0,0] = 546.67
A[0,1] = -137.82
A[0,2] = 94.083
A[0,3] = -204.68
A[0,4] = 111.51

A[1,0] = -1745.6
A[1,1] = 627.22
A[1,2] = -1175.2
A[1,3] = 2445.7
A[1,4] = -1341.7

A[2,0] = 3928.8
A[2,1] = -1519.3
A[2,2] = 4015.8
A[2,3] = -8415.3
A[2,4] = 4642.1

A[3,0] = -4384.8
A[3,1] = 1748.7
A[3,2] = -5362.1
A[3,3] = 11257.0
A[3,4] = -6218.2

A[4,0] = 1842.3
A[4,1] = -765.53
A[4,2] = 2507.7
A[4,3] = -5210.7
A[4,4] = 2867.5

result = 0.0

for i=0,4 do begin
for j=0, 4 do begin

    result = result+( A[i,j]* x^i * y^j)

endfor
endfor

result = result * Omega_m

;-------------------------------------------------------------------------
; Bryan
;-------------------------------------------------------------------------

if keyword_set(bryan) then begin
    x = Omega_m-1
    result = (18*!pi^2)+(82*x)-(39*x^2)  ; Scaling factor on critical density to match virial radius
endif

return,result
end
