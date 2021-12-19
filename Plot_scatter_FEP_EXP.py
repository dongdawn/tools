#! /usr/bin/python
# -*- coding: utf-8 -*-
# @Date    : 2021-12-18 14:05:51
# @Author  : ZYZCJHWDD 
# @Link    : https://github.com/dongdawn
# @Version : v1
import os
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
def rsquared(x, y, degree=1):
    results = {}
    coeffs = np.polyfit(x, y, degree)
    results['polynomial'] = coeffs.tolist()
    p = np.poly1d(coeffs)
    yhat = p(x)
    ybar = np.sum(y) / len(y)
    ssreg = np.sum((yhat - ybar) ** 2)
    sstot = np.sum((y - ybar) ** 2)
    results['determination'] = ssreg / sstot
    return results['determination']

def rmse(a, b):
    cnt = 0
    for i in range(len(a)):
        cnt += (a[i] - b[i]) ** 2
    return np.sqrt(cnt / len(a))

def mae(a, b):
    cnt = 0
    for i in range(len(a)):
        cnt += np.abs(a[i] - b[i])
    return cnt / len(a)
font_path = '/home/dongdong/wdd/calibribold.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=25)
title_prop = font_manager.FontProperties(fname=font_path, size=21)
def corr_plot(x_mean, y_mean, x_std, y_std, title=None, xlabel=None, ylabel=None, outputfile="corr_plot.png"):
    '''
    Plot correlation
    x_mean : list or numpy.array , mean of x
    y_mean : list or numpy.array , mean of y
    x_std : list or numpy.array , standard deviation of x
            if x_std is None, x_std will be set to zero list, whose length is equal to x_mean.
    y_std : list or numpy.array , standard deviation of y
    title : str
    xlabel : str
    ylabel : str
    outputfile : str , the filename of output.
    '''
    if x_std is None:
        x_std = [0 for ii in range(len(x_mean))]
    fig, ax = plt.subplots(figsize=(6, 6))
    xmin=min(x_mean)-0.1*(max(x_mean)-min(x_mean))
    ymin=min(y_mean)-0.1*(max(y_mean)-min(y_mean))
    xmax=max(x_mean)+0.1*(max(x_mean)-min(x_mean))
    ymax=max(y_mean)+0.1*(max(x_mean)-min(x_mean))
    r_min = min(xmin, ymin)
    r_max = max(xmax, ymax)
    diagline = np.linspace(r_min, r_max, 20)
    ax.fill_between(diagline, diagline - 1 * 2, diagline + 1 * 2, color='grey', alpha=0.2, zorder=0)
    ax.fill_between(diagline, diagline - 1 * 1, diagline + 1 * 1, color='grey', alpha=0.4, zorder=1)
    ax.plot(diagline, diagline, color='black', linewidth=3.)
    ax.errorbar(x_mean, y_mean, xerr=x_std, yerr=y_std, linestyle='None', marker='X',ms=12)
    ax.set_xlim(r_min, r_max)
    ax.set_ylim(r_min, r_max)
    if xlabel:
        #ax.set_xlabel(xlabel, fontsize=24)
        ax.set_xlabel(xlabel, fontproperties=font_prop)
    if ylabel:
        ax.set_ylabel(ylabel, fontproperties=font_prop)
    if title:
        ax.set_title(title, fontproperties=title_prop)
    ax.set_aspect('equal')
    fig.tight_layout()
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontproperties(font_prop)
    fig.savefig(outputfile, dpi=fig.dpi,bbox_inched='tight',pad_inches = 3)
    #plt.show()
    

data=np.loadtxt('fep_results.cs')

x_mean=data[:,0]
y_mean=data[:,1]
title = r"N : %d , RMSE : %.2f , R$^2$ : %.2f " % (len(x_mean), rmse(x_mean, y_mean), rsquared(x_mean, y_mean))
print(rsquared(x_mean, y_mean))
corr_plot(data[:,0], data[:,1],x_std=None, y_std=data[:,2], title=title, xlabel=r'$\Delta{G_{exp}}$(kcal/mol)', ylabel=r'$\Delta{G_{FEP}}$(kcal/mol)', outputfile="corr_plot.png")
