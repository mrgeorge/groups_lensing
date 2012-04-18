#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import fitsio
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt
import numpy

file="~/data/cosmos/code/group5_20110914.fits"
group=fitsio.read(file, ext=1)

# exclude flagged groups
group=group[group['FLAG_INCLUDE']==1]

msmr=group['ID_MMGG_SCALE'] != group['ID_MMGG_R200']
msbs=group['ID_MMGG_SCALE'] != group['ID_MLGG_SCALE']

# use helvetica and latex
plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20})
plt.rc('text', usetex=True)
plt.rc('axes',linewidth=1.5)

# start the plot with axes
plt.figure(1)
zmax=1.03
plt.xlim((0,zmax))
plt.ylim((41.2,44))
plt.xlabel(r'z',fontsize='medium')
plt.ylabel(r'log(L$_{\mathrm{X}}$) (erg/s)',fontsize='medium')

label1=r'$\mathrm{MMGG_{scale} \ne MMGG_{R200}}$'
label2=r'$\mathrm{MMGG_{scale} \ne BGG_{scale}}$'

plt.scatter(group['ZPHOT'],group['LX_APP'],c='gray',edgecolors='none',marker='o',s=20)
plt.scatter(group['ZPHOT'][msmr],group['LX_APP'][msmr],c='r',edgecolors='maroon',facecolors='r',marker='d',s=140,label=label1)
plt.scatter(group['ZPHOT'][msbs],group['LX_APP'][msbs],c='b',edgecolors='skyblue',facecolors='b',marker='s',s=70,label=label2)

plt.legend((label1, label2), loc='lower right', scatterpoints=1, markerscale=1, prop={'size':18},handletextpad=0).draw_frame(False)
ltext=plt.gca().get_legend().get_texts()
plt.setp(ltext[0],color='red')
plt.setp(ltext[1],color='blue')

# read in luminosity limits
groupLimFile="/Users/mgeorge/data/cosmos/auxfiles/lx_z_flm15p0.hst" # shallower limit for 96% of field, inner 52% goes deeper
zlim,llim=numpy.loadtxt(groupLimFile,unpack=True)
zlt1=zlim <= zmax
# 3rd order polynomial fit to measured curve
poly=numpy.polyfit(zlim,llim,3)
plt.plot(zlim,poly[0]*zlim**3 + poly[1]*zlim**2 + poly[2]*zlim + poly[3],'k--',linewidth=2)


plotDir="/Users/mgeorge/data/cosmos/groups_lensing/plots/"
plotFile=plotDir+"bad_groups.eps"
plt.savefig(plotFile)
