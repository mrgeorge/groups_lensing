 function nfw_sig, r, p, zl, $
                 r200 = r200, $
                 r180 = r180

;-------------------------------------------------------------------------
; BASED on Dave's : nfw_sigma.pro
; Wright & Brainerd
; 'Gravitational Lensing by NFW Haloes, 2000'
; Returns Sigma in units of h Msun pc^-2
;
; r200   -> 200 x critical density
; r180   -> 180 x mean density
; rvir   -> Delta_vir x critical density
;
; r is distance in MPC
; p[0] is r200 the virial (or 200) radius in MPC
; p[1] is c the concentration parameter
;-------------------------------------------------------------------------

defsysv,'!Omega_m', exists=exists
if NOT exists THEN define_cosmo

r_nfw  = p[0]
c      = p[1]

;-------------------------------------------------------------------------
; Calculate the critical density at z = zl
;-------------------------------------------------------------------------

rho_crit = critical_density(zl)

;-------------------------------------------------------------------------
; Overdensity
; (default definition is virial)
;-------------------------------------------------------------------------

overdensity = get_overdensity(zl)  ; Virial
if keyword_set(r200) then overdensity = get_overdensity(zl,/r200)
if keyword_set(r180) then overdensity = get_overdensity(zl,/r180)

;-------------------------------------------------------------------------
; NFW STUFF
;-------------------------------------------------------------------------

;CHARACTERISTIC OVERDENSITY
del_c = overdensity * (1.0/3.0)*c^3 /(alog(1.0+c)-c/(1.0+c))

;SCALE RADIUS
rs = (r_nfw*1.0)/(c*1.0)     ; rs = r_nfw/c

xx  = double(r/rs)           ; Parametrized as funct r/rs 
x   = 0.0

; FACTOR (cf paper)
; there is a factor of 2 here as opposed to DS
factor = 2*rs*del_c*rho_crit

;three cases to consider
ep = double(0.001)                          ;buffer , use linear interp in here

sig=dblarr(n_elements(r))
w1=where(xx lt (1.0-ep),n1)
w2=where(xx ge (1.0-ep) and xx le (1.0+ep),n2)
w3=where(xx gt (1.0+ep),n3)

if n1 gt 0 then begin
    x=xx[w1]
    s=(1 -2/sqrt(1 - x^2)*atanh(sqrt((1-x)/(1+x))))/(x^2-1)
    sig[w1]=s
endif

if n3 gt 0 then begin
    x=xx[w3]
    s=(1 -2/sqrt(x^2 -1)*atan(sqrt((x-1)/(1+x))))/(x^2-1)
    sig[w3]=s
endif

if n2 gt 0 then begin
    x=xx[w2]
    s=1/3d                     
    x1=1.0-ep
    x2=1.0+ep
    s1=(1 -2/sqrt(1 - x1^2)*atanh(sqrt((1-x1)/(1+x1))))/(x1^2-1)
    s2=(1 -2/sqrt(x2^2 -1)*atan(sqrt((x2-1)/(1+x2))))/(x2^2-1)
    s=(x-x1)*s2/(x2-x1) + (x-x2)*s1/(x1-x2)
    sig[w2]=s
endif

; In units of h Msun Mpc-2
sig = sig*factor

; Convert into units of h Msun pc^-2
sig = sig*1e-12

return,sig
end
