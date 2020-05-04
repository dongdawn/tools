#!/usr/bin/python
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
font_path = '/home/disk2/wdd/msm/calibribold.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=20)
lab_prop = font_manager.FontProperties(fname=font_path, size=15)

lagtime = 110        #fs
population=np.loadtxt('Populations.dat')
dataMSM=np.loadtxt('ProbMSM')
dataMD=np.loadtxt('ProbMD')
dataMDvar=np.loadtxt('var')
sortIndex=np.argsort(-population)
print population[sortIndex]
num = len(population)
numck = len(dataMSM[0])
x=np.arange(0,numck)*lagtime
fig = plt.figure(figsize=(6.6,4.8))
for i in range(9):
    sub = fig.add_subplot(3,3,i+1)
    MSM=dataMSM[sortIndex[i],:]
    MD=dataMD[sortIndex[i],:]
    var=dataMDvar[sortIndex[i],:]
    sub.errorbar(x, MD, yerr=var, fmt='o', label='MD',markerfacecolor='none')
    sub.plot(x, MSM, '-', label='MSM',lw=2.2)
    p = '{:.2e}'.format(population[sortIndex[i]])
    labP='P = ' +str(p)
    #plt.text(lagtime*0.3,1.02,labP,fontproperties=lab_prop)
    plt.xlim(-lagtime/2,lagtime*(numck-1)*1.11)
    plt.ylim(-0.1,1.3)
    sub.set_yticks(np.linspace(0,1,3))
    sub.set_xticks(np.linspace(0,lagtime*3,4))
    sub.set_xticklabels([])
    sub.set_yticklabels([])
    if i+1==1 or i+1==4 or i+1==7:
        sub.set_yticklabels([0.0,0.5,1.0])
    if i+1==7 or i+1==8 or i+1==9:
        sub.set_xticklabels([0,lagtime*1,lagtime*2,lagtime*3])
    if i+1==4:
        plt.ylabel('residence probability',fontproperties=font_prop) 
    if i+1==8:
        plt.xlabel('lag time (ns)',fontproperties=font_prop) 
    for label in (sub.get_xticklabels() + sub.get_yticklabels()):
        label.set_fontproperties(font_prop)
        label.set_fontsize(15)
plt.savefig('ckplot.png',dpi=600,bbox_inches='tight')
plt.show()
