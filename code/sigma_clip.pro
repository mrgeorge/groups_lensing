pro sigma_clip,array,amean,asigma,nsig=nsig,nIter=nIter,$
print=print,plot=plot,_extra=ex,pause=pause,index=index

;a simple sigma clipping algoritm

if n_params() eq 0 or not keyword_set(nsig) or $
	not keyword_set(nIter) then begin
	print, '-syntax sigma_clip,arr,mean,sigma,nsig=nsig,'
	print,'nIter=nIter,print=print,pause=pause,_extra=ex,plot=plot'
	return
endif

arr=array
wif=n_elements(arr)
index=lindgen(wif)
parts=(max(array)-min(array))/50
for i=0,nIter-1 do begin
if keyword_set(plot) then begin
	if i eq 0 then plothist,arr,_extra=ex,bin=parts
	if i gt 0 then plothist2,arr,/overplot,bin=parts,_extra=ex
endif

mom=moment(arr)
m=mom(0)
s=sqrt(mom(1))
clip=nsig*s
if keyword_set(print) then print,' mean',m,'   sigma',s,'  number',wif 
w=where(abs(arr-m) lt clip,wif)
if wif eq 0 then begin
	print,'stoped in sigma_clip.pro: nsig is too small. Everything clipped on iteration',i
        stop
	mom=moment(arr)
	amean=mom(0)
	asigma=0
	return
endif 

if n_elements(pause) gt 0 then wait,pause
arr=arr(w)
index=index(w)
endfor

mom=moment(arr)
amean=mom(0)
asigma=sqrt(mom(1))
if keyword_set(print) then print,' final mean',m,'   final sigma',s
return
end


