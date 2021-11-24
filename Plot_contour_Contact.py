import os
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.colors as colors
import matplotlib as ml

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

cmap = plt.get_cmap('jet')
new_cmap = truncate_colormap(cmap, 0.05, 0.97)

font_path = '/home/disk2/wdd/msm/calibribold.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=20)
leg_prop = font_manager.FontProperties(fname=font_path, size=15.5)
lab_prop = font_manager.FontProperties(fname=font_path, size=20)
name = "contact_all.cs"
data = np.loadtxt(name)
H=np.zeros((67,67))
for line in data:
    H[line[0]][line[1]]=line[2]
print H
HH = (H+H.T)/2
#Hn = HH - np.amin(HH)
fig = plt.figure(figsize=(5.1, 4))
ax = fig.add_subplot(111)
#cmap = plt.cm.get_cmap("jet_r")
#cmap.set_over("white")
#CS = plt.contourf(range(1,68,1),range(1,68,1),HH,levels = np.linspace(0,100,11),cmap=cmap,extend="max")
CS = plt.contourf(range(1,68,1),range(1,68,1),HH,levels = np.linspace(0,90,11),cmap=new_cmap)
cbar = plt.colorbar(CS)
#plt.xlabel(r'hydrophobic ($\mathregular{\AA^2}$)',fontproperties=font_prop)
#plt.ylabel(r'hydrophilic ($\mathregular{\AA^2}$)',fontproperties=font_prop)
plt.xlabel(r'residue index',fontproperties=font_prop)
plt.ylabel(r'residue index',fontproperties=font_prop)
#leg=plt.legend(loc=1, labelspacing=0.1, prop=leg_prop, scatterpoints=1, markerscale=1, numpoints=1,handlelength=1.5)
#leg.get_frame().set_linewidth(0.0)
#leg.get_frame().set_alpha(0.1)
for label in (ax.get_xticklabels() + ax.get_yticklabels()+cbar.ax.get_yticklabels()):
    label.set_fontproperties(font_prop)
for label in cbar.ax.get_yticklabels():
    label.set_fontproperties(leg_prop)
plt.savefig('contact_all.png',dpi=600,bbox_inches='tight')
    #label.set_fontsize(16)
#plt.imshow(HH, interpolation='nearest', origin='low',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],levels = np.linspace(1,6,20))
#plt.show()
