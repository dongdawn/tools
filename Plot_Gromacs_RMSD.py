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
import pandas as pd
import seaborn as sns
sns.set(color_codes=True)
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

def parse_arg():
    parser = argparse.ArgumentParser(description='This is a program to plot rmsd as time')
    parser.add_argument('-infile', dest='infile', help="input dat, this file is from dssp program,do_dssp -s %s.tpr -f %s_fit.xtc -o dssp_%s.xpm -sc scount_%s.xvg -ssdump ssdump_%s.dat -tu ns ", default='ssdump.dat')
    parser.add_argument('-outfile', dest='outfile', help="output file", required=True)
    arg = parser.parse_args()
    return arg.infile, arg.outfile

def plot_rmsd(filename,figout):
    font_path = '/home/disk2/wdd/msm/calibribold.ttf'
    font_prop = font_manager.FontProperties(fname=font_path, size=24)
    leg_prop = font_manager.FontProperties(fname=font_path, size=17)
    os.system("sed -i 's/^@/#/g' %s " %filename)
    fig, ax = plt.subplots(figsize=(6,4.7))
    data = np.loadtxt(filename)
    time = data[:,0]
    rmsd = data[:,1]
    ax.plot(time,rmsd)
    ax.set_ylabel('rmsd (nm)',fontproperties=font_prop)
    ax.set_xlabel('time (ns)',fontproperties=font_prop)
    plt.ylim(0,np.max(rmsd)*1.2)
    plt.xlim(0,np.max(time))
    #leg=plt.legend(loc=1, labelspacing=0.1, prop=leg_prop, scatterpoints=1, markerscale=1, numpoints=1,handlelength=1.5)
    #leg.get_frame().set_linewidth(0.0)
    #leg.get_frame().set_alpha(0.1)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontproperties(font_prop)
        label.set_fontsize(16)
    plt.savefig(figout,dpi=600,bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    infile, outfile= parse_arg()
    plot_rmsd(infile,outfile)

