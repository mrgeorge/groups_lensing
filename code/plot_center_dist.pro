pro plot_center_dist, groupInFile, outPathName

; add trailing slash to output path if missing
if(strmid(outPathName,strlen(outPathName)-1,1) NE '/') then outDir=outPathName+'/' else outDir=outPathName

group=mrdfits(groupInFile,1)
group=group[where(group.flag_include EQ 1,ngroups)]

;centers=['xray','mlgg_scale','mlgg_r200','mmgg_scale','mmgg_r200','cl','cm']
;center_str=textoidl(['X-ray','L_{max}(R_{scale})','L_{max}(R_{200})','M_{max}(R_{scale})','M_{max}(R_{200})','CL','CM'])
centers=['mmgg_scale','mmgg_r200','mlgg_scale','mlgg_r200','cn','cl','cm','xray']
center_str=textoidl(['MMGG_{scale}','MMGG_{R200}','MLGG_{scale}','MLGG_{R200}','CN','CL','CM','X-ray'])

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
label=replicate('',30)

for i=0,ncenters-1 do begin
   get_center_coords, group, centers[i],ra1,dec1,sel1
   for j=i,ncenters-1 do begin

      if(j EQ i) then begin
         plot,[0,1],[0,1],xstyle=4,ystyle=4,position=[pleft+i*pwidth,pbott+i*pheight,pleft+(i+1)*pwidth,pbott+(i+1)*pheight],thick=5
         continue
      endif else begin
;         compare_ctr,centers[i],centers[j],x_as,y_as,x_mpc,y_mpc
         get_center_coords, group, centers[j],ra2,dec2,sel2

         d_as=distance(ra1,dec1,ra2,dec2)*3600.
         d_mpc=d_as*(group.lensing_r200_mpc/group.lensing_r200_as)
         plothist,d_as,x_as,y_as,/noplot,bin=10
         plothist,d_mpc,x_mpc,y_mpc,/noplot,bin=0.05

         if(i EQ 0) then ytickname=label else ytickname=nolabel
         xtickname=nolabel
         plot,[0,150],[0,1],/nodata,position=[pleft+i*pwidth,pbott+j*pheight,pleft+(i+1)*pwidth,pbott+(j+1)*pheight],xrange=[0,150],yrange=[0,1],xtickname=xtickname,ytickname=ytickname,color=0,xstyle=1
         oplot,x_as,y_as/float(ngroups),psym=10,thick=5,color=!red
         
         if(j EQ ncenters-1) then axis,xaxis=1,xtickname=label

         if(i EQ 0) then xtickname=label else xtickname=nolabel
         ytickname=nolabel
         plot,[0,0.8],[0,1],/nodata,position=[pleft+j*pwidth,pbott+i*pheight,pleft+(j+1)*pwidth,pbott+(i+1)*pheight],xrange=[0,0.8],yrange=[0,1],xtickname=xtickname,ytickname=ytickname,color=0,xstyle=1
         oplot,x_mpc,y_mpc/float(ngroups),psym=10,thick=5,color=!blue
         
         if(j EQ ncenters-1) then axis,yaxis=1,ytickname=label

      endelse
   endfor
   
   xyouts,0.02,pbott+pheight/2.+i*pheight,center_str[i],/normal,alignment=0.5,orientation=90

   xyouts,pleft+pwidth/2.+i*pwidth,0.01,center_str[i],/normal,alignment=0.5
endfor

xyouts,pleft/2.+0.01,pbott+ncenters*pheight/2.,textoidl('% of groups'),/normal,alignment=0.5,orientation=90
xyouts,pleft+ncenters*pwidth/2.,pbott/2.,textoidl('Radial offset (Mpc)'),/normal,alignment=0.5
xyouts,pleft+ncenters*pwidth/2.,1-ptop/2.,textoidl('Radial offset (arcsec)'),/normal,alignment=0.5

device,/close


;------------------------------
; Make individual offset plots
;------------------------------
centers=['xray','mlgg_scale','mlgg_r200','mmgg_r200','cl','cm']
refcen='mmgg_scale'
get_center_coords, group, refcen,ra1,dec1,sel1

!p.multi=0

for i=0,n_elements(centers)-1 do begin
   get_center_coords, group, centers[i],ra2,dec2,sel2
   d_as=distance(ra1,dec1,ra2,dec2)*3600.
   d_kpc=d_as*(group.lensing_r200_mpc*1000./group.lensing_r200_as)

   device,filename=outDir+'dist_'+refcen+'_'+centers[i]+'.eps',/encapsul,/color,/helvetica,/inches,xsize=6,ysize=4

   plothist,d_kpc,bin=50,xrange=[0,600],xstyle=1,yrange=[0,60],ystyle=1,xtitle='Radial Offset (kpc)',ytitle=textoidl('N_{groups}')

   device,/close
endfor

end
