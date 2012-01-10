pro plot_sample_models

fitType = [$
          1,$                   ; 0  M0    : baryonic mass
          1,$                   ; 1  R_vir : NFW virial mass
          1,$                   ; 2  C     : NFW concentration
          0,$                   ; 3  alpha : fraction
          0,$                   ; 4  bias
          0,$                   ; 5  m_sigma
          1]                    ; 6  offset

fidM0=10.8
fidMh=13.4
fidConc=4.
fidRoff=0.

p_mean1=[fidM0,fidMh,fidConc,fidRoff]
p_mean2=[10.0,fidMh,fidConc,fidRoff]
p_mean3=[11.6,fidMh,fidConc,fidRoff]
p_mean4=[fidM0,13.1,fidConc,fidRoff]
p_mean5=[fidM0,13.7,fidConc,fidRoff]
p_mean6=[fidM0,fidMh,2.,fidRoff]
p_mean7=[fidM0,fidMh,10.,fidRoff]
p_mean8=[fidM0,fidMh,fidConc,0.05]
p_mean9=[fidM0,fidMh,fidConc,0.1]

zl=0.5
msun=11.3
x_mpc=0.01*10.^(findgen(100)/99*2.5)
cen_type='rhotis'
off_type='max3d'

get_ds_model,fitType,p_mean1,zl,msun,x_mpc,/use_m200,cen_type=cen_type,off_type=off_type,ps_term=fidCen,nfw_term=fidNFW
get_ds_model,fitType,p_mean2,zl,msun,x_mpc,/use_m200,cen_type=cen_type,off_type=off_type,ps_term=M0Lo
get_ds_model,fitType,p_mean3,zl,msun,x_mpc,/use_m200,cen_type=cen_type,off_type=off_type,ps_term=M0Hi
get_ds_model,fitType,p_mean4,zl,msun,x_mpc,/use_m200,cen_type=cen_type,off_type=off_type,nfw_term=MhLo
get_ds_model,fitType,p_mean5,zl,msun,x_mpc,/use_m200,cen_type=cen_type,off_type=off_type,nfw_term=MhHi
get_ds_model,fitType,p_mean6,zl,msun,x_mpc,/use_m200,cen_type=cen_type,off_type=off_type,nfw_term=cLo
get_ds_model,fitType,p_mean7,zl,msun,x_mpc,/use_m200,cen_type=cen_type,off_type=off_type,nfw_term=cHi
get_ds_model,fitType,p_mean8,zl,msun,x_mpc,/use_m200,cen_type=cen_type,off_type=off_type,nfw_term=offLo
get_ds_model,fitType,p_mean9,zl,msun,x_mpc,/use_m200,cen_type=cen_type,off_type=off_type,nfw_term=offHi

!p.multi=0
!p.font=0
!p.thick=8
!x.thick=5
!y.thick=5
!p.charsize=1.2
!p.charthick=1.2

set_plot,'ps'
simpctable
sun=sunsymbol()

xr=[0.03,3]
yr=[0,200.]

device,filename='sampleM0.eps',/color,/helvetica,/encapsul

plot,/nodata,xr,yr,xst=1,yst=1,/xlog,xtitle='Physical transverse radius R (Mpc)',ytitle=textoidl('\Delta\Sigma(R) (M')+sun+textoidl(' pc^{-2})')
oplot,x_mpc,fidNFW
oplot,x_mpc,fidCen,color=!green,linestyle=2
oplot,x_mpc,M0Lo,color=!blue,linestyle=2
oplot,x_mpc,M0Hi,color=!red,linestyle=2

device,/close
device,filename='sampleMh.eps',/color,/helvetica,/encapsul

plot,/nodata,xr,yr,xst=1,yst=1,/xlog,xtitle='Physical transverse radius R (Mpc)',ytitle=textoidl('\Delta\Sigma(R) (M')+sun+textoidl(' pc^{-2})')
oplot,x_mpc,fidCen,linestyle=2,color=!green
oplot,x_mpc,fidNFW
oplot,x_mpc,MhLo,color=!blue
oplot,x_mpc,MhHi,color=!red

device,/close
device,filename='sampleC.eps',/color,/helvetica,/encapsul

plot,/nodata,xr,yr,xst=1,yst=1,/xlog,xtitle='Physical transverse radius R (Mpc)',ytitle=textoidl('\Delta\Sigma(R) (M')+sun+textoidl(' pc^{-2})')
oplot,x_mpc,fidCen,linestyle=2,color=!green
oplot,x_mpc,fidNFW
oplot,x_mpc,cLo,color=!blue
oplot,x_mpc,cHi,color=!red

device,/close
device,filename='sampleRoff.eps',/color,/helvetica,/encapsul

plot,/nodata,xr,yr,xst=1,yst=1,/xlog,xtitle='Physical transverse radius R (Mpc)',ytitle=textoidl('\Delta\Sigma(R) (M')+sun+textoidl(' pc^{-2})')
oplot,x_mpc,fidCen,linestyle=2,color=!green
oplot,x_mpc,fidNFW
oplot,x_mpc,offLo,color=!blue
oplot,x_mpc,offHi,color=!red

device,/close

end
