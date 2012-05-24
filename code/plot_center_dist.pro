pro pad_hist,x=x,y=y
; add halfbins to both ends of a histogram
; x,y are e.g. the outputs from plothist
; returned with values appended

; note, this connects the histogram with the y-axis on the upper left side
; of the plot and the x-axis on the bottom right
; unless the left-most point is still >0, then it connects it down to
; the x-axis.

xstep=x[1]-x[0]
x=[x[0]-0.5*xstep,x,replicate(x[n_elements(x)-1]+0.5*xstep,2)]
y=[y[0],y,y[n_elements(y)-1],0]
if(min(x) GT 0) then begin
   x=[x[0],x]
   y=[0,y]
endif

end

pro plot_center_dist, groupInFile, outPathName


if(n_elements(groupInFile) EQ 0) then groupInFile="~/data/cosmos/code/group5_cenunc_20110914_run20120518.fits"
if(n_elements(outPathName) EQ 0) then outPathName="~/data/cosmos/groups_lensing/plots/"

; add trailing slash to output path if missing
if(strmid(outPathName,strlen(outPathName)-1,1) NE '/') then outDir=outPathName+'/' else outDir=outPathName

group=mrdfits(groupInFile,1)
group=group[where(group.flag_include EQ 1,ngroups)]

;centers=['xray','mlgg_scale','mlgg_r200','mmgg_scale','mmgg_r200','cl','cm']
;center_str=textoidl(['X-ray','L_{max}(R_{scale})','L_{max}(R_{200})','M_{max}(R_{scale})','M_{max}(R_{200})','CL','CM'])
centers=['mmgg_scale','mmgg_r200','bgg_scale','bgg_r200','cn','cm','cf','xray']
center_str=textoidl(['MMGG_{scale}','MMGG_{R200}','BGG_{scale}','BGG_{R200}','CN','CM','CF','X-ray'])

ncenters=n_elements(centers)

!p.multi=[0,ncenters,ncenters]
!p.font=0
!p.thick=3
!x.thick=3
!y.thick=3
!p.charsize=1.2

set_plot,'ps'
simpctable
device,filename=outDir+'center_dist.eps',/helvetica,/color,xsize=10,ysize=10,/inches,/encapsul

pleft=0.08
pright=0.03
pbott=0.08
ptop=0.05
pwidth=(1.-pleft-pright)/ncenters
pheight=(1.-pbott-ptop)/ncenters

nolabel=replicate(' ',30)
ytickv=[0,0.2,0.4,0.6,0.8]
ytickname_frac=['0','0.2','0.4','0.6','0.8']
yticks=n_elements(ytickv)-1
yminor=2
xtickv_mpc=[0,0.2,0.4,0.6]
xtickname_mpc=['0','0.2','0.4','0.6'];string(xtickv_mpc,format='(F3.1)')
xticks_mpc=n_elements(xtickv_mpc)-1
xminor_mpc=4
xtickv_as=[0,50,100]
xtickname_as=strcompress(string(xtickv_as,format='(I)'),/remove)
xticks_as=n_elements(xtickv_as)-1
xminor_as=5

xrange_as=[0,150]
xrange_mpc=[0,0.75]
yrange=[0,1.0]

axischarsize=2.

for i=0,ncenters-1 do begin
   get_center_coords, group, centers[i],ra1,dec1,sel1
   for j=i,ncenters-1 do begin

      if(j EQ i) then begin
         ; plot uncertainty distribution for centroids, else blank
         if(centers[i] EQ 'xray') then begin
            plothist,group.pos_err_ellipse*3600.*group.lensing_r200_mpc/group.lensing_r200_as,x_mpc,y_mpc,bin=0.05,/noplot
            plothist,group.pos_err_ellipse*3600.*group.lensing_r200_mpc/group.lensing_r200_as,x_mpc,y_mpc,bin=0.05,peak=max(y_mpc)/float(ngroups),xrange=xrange_mpc,yrange=yrange,position=[pleft+i*pwidth,pbott+(ncenters-1-j)*pheight,pleft+(i+1)*pwidth,pbott+(ncenters-j)*pheight],xtickname=xtickname_mpc,xtickv=xtickv_mpc,xticks=xticks_mpc,xminor=xminor_mpc,ytickname=nolabel,yticks=yticks,ytickv=ytickv,yminor=yminor,xstyle=1,ystyle=1,charsize=axischarsize,/fill,color=!purple,fcolor=!purple
            axis,yaxis=1,yrange=yrange,ystyle=1,ytickname=ytickname_frac,ytickv=ytickv,yticks=yticks,yminor=yminor,charsize=axischarsize
         endif else if(centers[i] EQ 'cn') then begin
            plothist,group.pos_err_cn*3600.*group.lensing_r200_mpc/group.lensing_r200_as,x_mpc,y_mpc,bin=0.05,/noplot
            plothist,group.pos_err_cn*3600.*group.lensing_r200_mpc/group.lensing_r200_as,x_mpc,y_mpc,bin=0.05,peak=max(y_mpc)/float(ngroups),xrange=xrange_mpc,yrange=yrange,position=[pleft+i*pwidth,pbott+(ncenters-1-j)*pheight,pleft+(i+1)*pwidth,pbott+(ncenters-j)*pheight],xtickname=nolabel,xtickv=xtickv_mpc,xticks=xticks_mpc,xminor=xminor_mpc,ytickname=nolabel,yticks=yticks,ytickv=ytickv,yminor=yminor,xstyle=1,ystyle=1,charsize=axischarsize,/fill,color=!purple,fcolor=!purple
         endif else if(centers[i] EQ 'cl' OR centers[i] EQ 'cf') then begin
            plothist,group.pos_err_cl*3600.*group.lensing_r200_mpc/group.lensing_r200_as,x_mpc,y_mpc,bin=0.05,/noplot
            plothist,group.pos_err_cl*3600.*group.lensing_r200_mpc/group.lensing_r200_as,x_mpc,y_mpc,bin=0.05,peak=max(y_mpc)/float(ngroups),xrange=xrange_mpc,yrange=yrange,position=[pleft+i*pwidth,pbott+(ncenters-1-j)*pheight,pleft+(i+1)*pwidth,pbott+(ncenters-j)*pheight],xtickname=nolabel,xtickv=xtickv_mpc,xticks=xticks_mpc,xminor=xminor_mpc,ytickname=nolabel,yticks=yticks,ytickv=ytickv,yminor=yminor,xstyle=1,ystyle=1,charsize=axischarsize,/fill,color=!purple,fcolor=!purple
         endif else if(centers[i] EQ 'cm') then begin
            plothist,group.pos_err_cm*3600.*group.lensing_r200_mpc/group.lensing_r200_as,x_mpc,y_mpc,bin=0.05,/noplot
            plothist,group.pos_err_cm*3600.*group.lensing_r200_mpc/group.lensing_r200_as,x_mpc,y_mpc,bin=0.05,peak=max(y_mpc)/float(ngroups),xrange=xrange_mpc,yrange=yrange,position=[pleft+i*pwidth,pbott+(ncenters-1-j)*pheight,pleft+(i+1)*pwidth,pbott+(ncenters-j)*pheight],xtickname=nolabel,xtickv=xtickv_mpc,xticks=xticks_mpc,xminor=xminor_mpc,ytickname=nolabel,yticks=yticks,ytickv=ytickv,yminor=yminor,xstyle=1,ystyle=1,charsize=axischarsize,/fill,color=!purple,fcolor=!purple
         endif else continue
      endif else begin
;         compare_ctr,centers[i],centers[j],x_as,y_as,x_mpc,y_mpc
         get_center_coords, group, centers[j],ra2,dec2,sel2

         d_as=distance(ra1,dec1,ra2,dec2)*3600.
         d_mpc=d_as*(group.lensing_r200_mpc/group.lensing_r200_as)

         plothist,d_as,x_as,y_as,/noplot,bin=10
         plothist,d_mpc,x_mpc,y_mpc,/noplot,bin=0.05

         pad_hist,x=x_as,y=y_as
         pad_hist,x=x_mpc,y=y_mpc

         ; plot the distribution of offsets in physical radii (bottom left)
         if(i EQ 0) then ytickname=ytickname_frac else ytickname=nolabel
         if(j EQ ncenters-1) then xtickname=xtickname_mpc else xtickname=nolabel

         plot,xrange_mpc,yrange,/nodata,position=[pleft+i*pwidth,pbott+(ncenters-1-j)*pheight,pleft+(i+1)*pwidth,pbott+(ncenters-j)*pheight],xtickname=xtickname,xtickv=xtickv_mpc,xticks=xticks_mpc,xminor=xminor_mpc,ytickname=ytickname,yticks=yticks,ytickv=ytickv,yminor=yminor,color=0,xstyle=1,ystyle=1,charsize=axischarsize
         oplot,x_mpc,y_mpc/float(ngroups),psym=10,thick=5,color=!blue

         ; plot the distribution of offsets in arcsec (upper right)
         xtickname=nolabel
         ytickname=nolabel
         plot,xrange_as,yrange,/nodata,position=[pleft+j*pwidth,pbott+(ncenters-1-i)*pheight,pleft+(j+1)*pwidth,pbott+(ncenters-i)*pheight],xtickname=xtickname,xticks=xticks_as,xtickv=xtickv_as,xminor=xminor_as,ytickname=ytickname,yticks=yticks,ytickv=ytickv,yminor=yminor,color=0,xstyle=1,ystyle=1,charsize=axischarsize
         oplot,x_as,y_as/float(ngroups),psym=10,thick=5,color=!red
         if(i EQ 0) then axis,xaxis=1,xrange=xrange_as,xstyle=1,xtickname=xtickname_as,xticks=xticks_as,xtickv=xtickv_as,xminor=xminor_as,charsize=axischarsize
         if(j EQ ncenters-1) then axis,yaxis=1,yrange=yrange,ystyle=1,ytickname=ytickname_frac,ytickv=ytickv,yticks=yticks,yminor=yminor,charsize=axischarsize
      endelse
   endfor
   
   ; label rows on left side
   xyouts,0.02,pbott+pheight/2.+(ncenters-1-i)*pheight,center_str[i],/normal,alignment=0.5,orientation=90
   ; label columns on bottom
   xyouts,pleft+pwidth/2.+i*pwidth,0.01,center_str[i],/normal,alignment=0.5
endfor

xyouts,pleft/2.+0.005,pbott+ncenters*pheight/2.,textoidl('Fraction of groups'),/normal,alignment=0.5,orientation=90
xyouts,pleft+ncenters*pwidth/2.,pbott/2.-0.005,textoidl('Physical transverse separation (h_{72}^{-1} Mpc)'),/normal,alignment=0.5
xyouts,pleft+ncenters*pwidth/2.,1-ptop/2.,textoidl('Angular separation (arcsec)'),/normal,alignment=0.5

device,/close


;------------------------------
; Make individual offset plots
;------------------------------
;centers=['xray','mlgg_scale','mlgg_r200','mmgg_r200','cl','cm']
;refcen='mmgg_scale'
;get_center_coords, group, refcen,ra1,dec1,sel1
;
;!p.multi=0
;
;for i=0,n_elements(centers)-1 do begin
;   get_center_coords, group, centers[i],ra2,dec2,sel2
;   d_as=distance(ra1,dec1,ra2,dec2)*3600.
;   d_kpc=d_as*(group.lensing_r200_mpc*1000./group.lensing_r200_as)
;
;   device,filename=outDir+'dist_'+refcen+'_'+centers[i]+'.eps',/encapsul,/color,/helvetica,/inches,xsize=6,ysize=4
;
;   plothist,d_kpc,bin=50,xrange=[0,600],xstyle=1,yrange=[0,60],ystyle=1,xtitle='Radial Offset (kpc)',ytitle=textoidl('N_{groups}')
;
;   device,/close
;endfor

end
