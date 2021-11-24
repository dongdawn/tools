#!/usr/bin/python
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from scipy.stats.stats import pearsonr
font_path = './calibribold.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=18)
wfexp=open('experiment_CA.cs','r')
wfmd=open('reweight_results.cs','r')
dataexp=wfexp.readlines()
datamd=wfmd.readlines()
expc=[]
mdc=[]
for i in range(67):
    linee=dataexp[i].strip()
    linem=datamd[i].strip()
    if linee!='****' and linem!='****':
        expc.append(float(linee))
        mdc.append(float(linem))
pearson,p_v=pearsonr(expc,mdc)
print round(pearson,3)
fig = plt.figure(figsize=(6,5.2))
sub = fig.add_subplot(1,1,1)
la='R = '+str(round(pearson,3))
sub.scatter(expc,mdc,c='dodgerblue',linewidth=0,s=35)
sub.plot([0,100],[0,100],'--',c='silver',alpha=0.6)
plt.text(41,67,la,fontproperties=font_prop)
#leg=plt.legend(loc='upper left', labelspacing=0.1, prop=font_prop, scatterpoints=1, markerscale=1, numpoints=1)
#leg.get_frame().set_linewidth(0.0)
#leg.get_frame().set_alpha(0.1)
plt.xlim(40,70)
plt.ylim(40,70)
plt.grid(alpha=0.5)
#sub.set_yticks(np.linspace(0,1,3))
#sub.set_xticks(np.linspace(0,1,3))
plt.ylabel('MD (ppm)',fontproperties=font_prop)
plt.xlabel('Expt. (ppm)',fontproperties=font_prop)
for label in (sub.get_xticklabels() + sub.get_yticklabels()):
    label.set_fontproperties(font_prop)
plt.savefig('CA.png',dpi=600,bbox_inches='tight')
plt.show()

