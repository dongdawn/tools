#! /usr/bin/python
# -*- coding: utf-8 -*-
# @Date    : 2018-01-30 14:05:51
# @Author  : WDD 
# @Link    : https://github.com/dongdawn
# @Version : v1
import sys
import os
import numpy as np
import argparse
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

def parse_arg():
    parser = argparse.ArgumentParser(description='This is a program to analysis the secondary structures')
    parser.add_argument('-infile', dest='infile', help="input dat, this file is from dssp program,do_dssp -s %s.tpr -f %s_fit.xtc -o dssp_%s.xpm -sc scount_%s.xvg -ssdump ssdump_%s.dat -tu ns ", default='ssdump.dat')
    parser.add_argument('-outfile', dest='outfile', help="output file", required=True)
    parser.add_argument('-plot_flag', dest='plot_flag', help="if plot",default='True')
    parser.add_argument('-figname', dest='figname', help="the output file name of plot")
    arg = parser.parse_args()
    return arg.infile, arg.outfile, arg.plot_flag, arg.figname

def statis(fn,outname):
    with open(fn,"r") as fin, open(outname,'w') as wf:
        data = fin.readline()
        datatmp = []
        for k in fin:
            datatmp.append(list(k.strip()))
        mat = np.array(datatmp)
        all = set(mat.flatten().tolist())
        wf.write("#COLUMN\t")
        for char in all:
            wf.write("%6s\t" %char)
        wf.write("\n")
        for i in range(mat.shape[1]):
            col = mat[:,i].tolist()
            wf.write("%6d\t" %(i+1))
            for char in all:
                wf.write("%6d\t" %col.count(char)) 
            wf.write("\n")

def plot_bar(statout,figout):
    font_path = '/home/disk2/wdd/msm/calibribold.ttf'
    font_prop = font_manager.FontProperties(fname=font_path, size=24)
    leg_prop = font_manager.FontProperties(fname=font_path, size=17)
    
    with open(statout,'r') as sto:
        firstLine = sto.readline()
    ss=firstLine.split()
    #print ss.index('H')
    width=0.4
    #fig = plt.figure(figsize=(6.5,5.5))
    fig, ax = plt.subplots(figsize=(6,4.7))
    data = np.loadtxt(statout)
    resid = data[:,0]
    num_frames=np.sum(data[0,1:])
    print num_frames
    helix = data[:,ss.index('H')]/num_frames
    if 'E' in ss:
        sheet = data[:,ss.index('E')]/num_frames
    else:
        sheet = [0]*len(data)
    coil = data[:,ss.index('~')]/num_frames
    turn = data[:,ss.index('T')]/num_frames
    
    bar_h = ax.bar(resid-width,helix,width, color='m',label=r'$\mathregular{\alpha}$-helix',edgecolor='m')
    bar_s = ax.bar(resid,sheet,width, color='gold',label=r'$\mathregular{\beta}$-sheet',edgecolor='gold')
    bar_c = ax.bar(resid+width,coil,width, color='silver',label='coil',edgecolor='silver')
    ax.set_ylabel('probability',fontproperties=font_prop)
    ax.set_xlabel('residue index',fontproperties=font_prop)
    plt.ylim(0,1.3)
    plt.xlim(0,len(data))
    leg=plt.legend(loc=1, labelspacing=0.1, prop=leg_prop, scatterpoints=1, markerscale=1, numpoints=1,handlelength=1.5)
    leg.get_frame().set_linewidth(0.0)
    leg.get_frame().set_alpha(0.1)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontproperties(font_prop)
        label.set_fontsize(16)
    plt.savefig(figout,dpi=600,bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    infile, outfile, plot_flag, figname = parse_arg()
    statis(infile,outfile)
    print outfile
    if plot_flag:
        plot_bar(outfile,figname)

