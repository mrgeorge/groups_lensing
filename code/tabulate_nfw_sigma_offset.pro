function tabulate_nfw_sigma_offset,r, roff, p, zl,$
                        r200 = r200, $
                        r180 = r180

; Tabulate nfw_sigma_offset. Valid between 1e-4 and 4Mpc
; Could make smaller by using a smaller array
; Note: in this function input r is one elements (not an array)

res = 0.0

common common_nfw_sig_offset, reset_tab, r_tab, nfw_sig_offset_tab 

;!EXCEPT=2

; Define variables if not previously defined
s=size(r_tab)

if(s[1] eq 0) then begin ; Tabulate r in log space
   npts = 2000.0 ; high number just to begin with
   fac = 1.02    ; 
   reset_tab = 1
   
   r_tab    = dindgen(npts)
   r_tab[0] = 1e-4  ; minimum valid radius
   for j=1,n_elements(r_tab)-1 do begin
      r_tab[j] = r_tab[j-1]*fac     
   endfor

   ; Valid to 4Mpc
   sel=where(r_tab le 4)
   r_tab = r_tab(sel)
   npts = n_elements(sel)
   nfw_sig_offset_tab = dindgen(npts)
endif

; check limits on r
if(r gt max(r_tab)) then begin
   print,'Stopped in tabulate_nfw_sigma_offset.pro, input r in out of limits (>3Mpc)'
   stop
endif

if (reset_tab eq 1) then begin  ; Retabulate
   for i=0,n_elements(r_tab)-1 do begin
      nfw_sig_offset_tab[i]=nfw_sigma_offset(r_tab[i],roff,p,zl,r200=keyword_set(r200),r180=keyword_set(r180))
   endfor
   reset_tab = 0
endif

; Calculate using tabulated values
sel=where(r_tab gt r)
index = sel[0]

if(index eq -1) then begin
   print,'Stopped in tabulate_nfw_sigma_offset.pro, invalid selection'
   stop
endif
 
; if close to this bin, return this value
if(abs(r-r_tab[index]) lt 1e-5) then begin
   res = nfw_sig_offset_tab[index]
endif else begin ; otherwise, interpolate between index and index-1, linear interpolation in log space is more precise
   if (r le min(r_tab)) then begin
      res = nfw_sig_offset_tab[0]
   endif else begin
       ; problems with this because DS goes neg
      ;a   = (alog10(nfw_sig_offset_tab[index])-alog10(nfw_sig_offset_tab[index-1]))/(alog10(r_tab[index])-alog10(r_tab[index-1]))
      ;b   = alog10(nfw_sig_offset_tab[index]) - (a*alog10(r_tab[index])) 
      ;res = 10.0^((a*alog10(r))+b)

      a   = ((nfw_sig_offset_tab[index])-(nfw_sig_offset_tab[index-1]))/((r_tab[index])-(r_tab[index-1]))
      b   = (nfw_sig_offset_tab[index]) - (a*(r_tab[index])) 
      res = (a*r)+b

      ;a=CHECK_MATH()
      ;if (a gt 0) then stop
   endelse
endelse
return, res
end
