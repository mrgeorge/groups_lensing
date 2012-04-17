pro plot_cartoon_models

; plot toy models for cartoon miscentering figure

plotDir="~/data/cosmos/groups_lensing/plots/"

fitType = [$
          0,$                   ; 0  M0    : baryonic mass
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

p_mean1=[fidMh,fidConc,fidRoff]
p_mean9=[fidMh,fidConc,0.05]

zl=0.5
msun=11.3
x_mpc=0.01*10.^(findgen(100)/99*2.5)
off_type='max3d'

get_ds_model,fitType,p_mean1,zl,msun,x_mpc,/use_m200,off_type=off_type,nfw_term=fidNFW
get_ds_model,fitType,p_mean9,zl,msun,x_mpc,/use_m200,off_type=off_type,nfw_term=offHi

!p.multi=0
!p.font=0
!p.thick=10
!x.thick=7
!y.thick=7
!p.charsize=2.5
!p.charthick=1.2

set_plot,'ps'
simpctable
sun=sunsymbol()
blank=replicate(' ',30)

xr=[0.03,3]
yr=[0,160.]

device,filename=plotDir+'cartoon_cen.eps',/color,/helvetica,/encapsul
plot,/nodata,xr,yr,xst=1+8,yst=1+8,xticks=1,yticks=1,xtickname=blank,ytickname=blank,/xlog,xtitle='log(R)',ytitle=textoidl('\Delta\Sigma(R)'),xmarg=5,ymarg=4
oplot,x_mpc,fidNFW,color=!blue
device,/close

device,filename=plotDir+'cartoon_off.eps',/color,/helvetica,/encapsul
plot,/nodata,xr,yr,xst=1+8,yst=1+8,xticks=1,yticks=1,xtickname=blank,ytickname=blank,/xlog,xtitle='log(R)',ytitle=textoidl('\Delta\Sigma(R)'),xmarg=5,ymarg=4
oplot,x_mpc,offHi,color=!purple
device,/close


stop
end
