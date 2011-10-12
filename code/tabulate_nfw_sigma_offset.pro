function tabulate_nfw_sigma_offset,r, roff, p, zl,$
                        r200 = r200, $
                        r180 = r180

; Tabulate nfw_sigma_offset. Valid between 1e-4 and 4Mpc
; Could make smaller by using a smaller array
; Note: in this function input r can be a scalar or array

res = 0.0

common common_nfw_sig_offset, reset_tab, r_tab, nfw_sig_offset_tab 

; Define variables if not previously defined
s=size(r_tab)

if(s[1] eq 0) then begin ; Tabulate r in log space
   reset_tab=1

   minRad=1.e-4 ; Mpc
   maxRad=4. ; Mpc
   npts=500 ; in log r
   step=(maxRad/minRad)^(1./npts)
   r_tab=minRad*step^indgen(npts)

   nfw_sig_offset_tab = dindgen(npts)
endif

; check limits on r
if(max(r) gt max(r_tab)) then begin
   print,'Stopped in tabulate_nfw_sigma_offset.pro, input r in out of limits (>3Mpc)'
   stop
endif

if (reset_tab eq 1) then begin  ; Retabulate
   nfw_sig_offset_tab=nfw_sigma_offset(r_tab,roff,p,zl,r200=keyword_set(r200),r180=keyword_set(r180))
   reset_tab = 0
endif

; Calculate using tabulated values
index=lonarr(n_elements(r))
res=fltarr(n_elements(r))
for ii=0,n_elements(r)-1 do index[ii]=min(where(r_tab GE r[ii]))

if(min(index) LT 0) then begin
   print,'Stopped in tabulate_nfw_sigma_offset.pro, invalid selection'
   stop
endif

; for r<minRad, fill with lowest tabulated ds value, else interpolate in log R
low=where(index EQ 0,nLow,complement=middle,ncomplement=nMiddle) ; low is where r<minRad, middle is where minRad<r<maxRad. r>maxRad has already been flagged
if(nLow GT 0) then res[low]=nfw_sig_offset_tab[0]
if(nMiddle GT 0) then res[middle]=interpol(nfw_sig_offset_tab,alog10(r_tab),alog10(r[middle]))

return, res
end
