#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import fitsio
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt

file="~/data/cosmos/code/group5_20110914.fits"
group=fitsio.read(file, ext=1)

# exclude flagged groups
group=group[group['FLAG_INCLUDE']==1]

msmr=group['ID_MMGG_SCALE'] != group['ID_MMGG_R200']
msbs=group['ID_MMGG_SCALE'] != group['ID_MLGG_SCALE']


plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':18})
plt.rc('text', usetex=True)
plt.figure(1)
plt.xlim((0,1.03))
plt.ylim((41.2,44))
plt.xlabel(r'z',fontsize='medium')
plt.ylabel(r'log(L$_{\mathrm{X}}$) (erg/s)',fontsize='medium')


label1=r'$\mathrm{MMGG_{scale} \ne MMGG_{R200}}$'
label2=r'$\mathrm{MMGG_{scale} \ne BGG_{scale}}$'

plt.scatter(group['ZPHOT'],group['LX_APP'],c='gray',marker='o',s=15,edgecolors='none')
plt.scatter(group['ZPHOT'][msmr],group['LX_APP'][msmr],edgecolors='red',marker='d',s=100,facecolors='none',label=label1)
plt.scatter(group['ZPHOT'][msbs],group['LX_APP'][msbs],edgecolors='blue',marker='s',s=80,facecolors='none',label=label2)

plt.legend((label1, label2), loc='lower right', scatterpoints=1, markerscale=1, prop={'size':16},handletextpad=0).draw_frame(False)
ltext=plt.gca().get_legend().get_texts()
plt.setp(ltext[0],color='red')
plt.setp(ltext[1],color='blue')


plotDir="/Users/mgeorge/data/cosmos/groups_lensing/plots/"
plotFile=plotDir+"bad_groups.eps"
plt.savefig(plotFile)
