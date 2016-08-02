#!/usr/bin/env python3

import sys, re, os.path
import pandas as pd
import matplotlib

import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pylab import *

params={'axes.linewidth' : .5}
rcParams.update(params)

# get the filename from the command line
filename = sys.argv[1]

f = open(filename);
fl = f.readlines();
f.close()


parline = -1

for idx, line in enumerate(fl):
    if re.match("^type.*",line) != None:
        parline = idx - 1;
        break;

# read in the csv file
if parline > 0:
    histdat = pd.read_csv(filename, nrows=parline-3, sep=";")
else:
    histdat = pd.read_csv(filename, sep=";")

# only take every tenth generation, otherwise too much data....
histdat = histdat[histdat["generation"] % 10 == 0]

# generate the figure

# initialize and specify size 
fig = plt.figure(figsize=(10,10))

num_rows = 5

# add first subplot
plt.subplot(num_rows,1,1)
plt.plot(histdat["generation"],histdat["meanz"],'b',
        histdat["generation"],histdat["epsilon"],'darkgreen'
        )
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'Elevation, $\bar{z}$')

# add second subplot
plt.subplot(num_rows,1,2)
plt.plot(histdat["generation"],histdat["meanm_m"],'b',
        histdat["generation"],histdat["meanm_g"],'r',
        histdat["generation"],histdat["meanm_e"],'magenta',
        linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'Maternal effects')
plt.legend((r'$m_{m}$',r'$m_{g}$',r'$m_{e}$'))

# add third subplot
plt.subplot(num_rows,1,3)
plt.plot(histdat["generation"],histdat["nsurv"],'peru',
        linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'Survivors')

# add fourth subplot
plt.subplot(num_rows,1,4)
plt.plot(
        histdat["generation"],histdat["meang"],'r',
        histdat["generation"],histdat["meanb"],'b',
        linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'Elevation, slope')
plt.legend((r'$g$',r'$b$'))

# add fourth subplot
plt.subplot(num_rows,1,5)
plt.plot(
        histdat["generation"],histdat["varg"],'peru',
        histdat["generation"],histdat["varm_m"],'darkgreen',
        histdat["generation"],histdat["varm_g"],'powderblue',
        histdat["generation"],histdat["varm_e"],'magenta',
        linewidth=1)
plt.ylabel(r'Variances')
plt.legend((r'$\sigma_{g}^{2}$',r'$\sigma_{m_{m}}^{2}$',r'$\sigma_{m_{g}}^{2}$',r'$\sigma_{m_{e}}^{2}$'))

graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf")
