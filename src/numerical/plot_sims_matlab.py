#!/usr/bin/env python3

import sys, re, os.path, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd

from pylab import *

params={'axes.linewidth' : .5}
rcParams.update(params)

# get the filename from the command line
filename = sys.argv[1]

# read in the csv file
histdat = pd.read_csv(filename,sep=";")

histdat = histdat[

# initialize and specify size 
fig = plt.figure(figsize=(10,18))

num_rows = 6

def to_log(row):
    return(pd.Series(math.log10(1 + row["time"])))

histdat["time"] = histdat.apply(to_log, axis=1)

# add first subplot
plt.subplot(num_rows,1,1)
plt.plot(histdat["time"],histdat["abar"],'g',linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on')
plt.ylabel(r'mean trait value, $\bar{x}$')
plt.legend((r'$\bar{a}$',r'$\bar{a}$'))

plt.subplot(num_rows,1,2)
plt.plot(histdat["time"],histdat["bbar"],'r',linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on')
plt.ylabel(r'mean reaction norm slope, $\bar{b}$')
plt.legend((r'$\bar{b}$',r'$\bar{b}$'))

plt.subplot(num_rows,1,3)
plt.plot(histdat["time"],histdat["mbar"],'b',linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on')
plt.ylabel(r'mean trait value, $\bar{m}$')
plt.legend((r'$\bar{m}$',r'$\bar{m}$'))

plt.subplot(num_rows,1,4)
plt.plot(histdat["time"],histdat["zbar"],'b',linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on')
plt.ylabel(r'mean phenotype, $\bar{z}$')

plt.subplot(num_rows,1,5)
plt.plot(histdat["time"],histdat["fitness"],'b',linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on')
plt.ylabel(r'mean fitness, $\bar{W}$')

plt.subplot(num_rows,1,6)
plt.plot(histdat["time"],histdat["epst"],'b',linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on')
plt.ylabel(r'proportion soft inheritance')
#plt.subplot(num_rows,1,2)
#plt.plot(histdat["time"],histdat["p00"],'g',histdat["time"],histdat["p01"],'y',histdat["time"],histdat["p10"],'r',histdat["time"],histdat["p11"],('#8900ff'))
#plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
#plt.ylabel(r'$\bar{p}_{ij}$')
#plt.legend((r'$p_{q_{1},e_{1}}$',r'$p_{q_{1},e_{2}}$',r'$p_{q_{2},e_{1}}$',r'$p_{q_{2},e_{2}}$'))
#plt.ylim(-0.05,1.05)
#
#plt.subplot(num_rows,1,3)
#plt.plot(histdat["time"],histdat["v00"],'g',histdat["time"],histdat["v01"],'y',histdat["time"],histdat["v10"],'r',histdat["time"],histdat["v11"],('#8900ff'))
#plt.ylabel(r'$\bar{v}_{ij}$')
#plt.legend((r'$v_{q_{1},e_{1}}$',r'$v_{q_{1},e_{2}}$',r'$v_{q_{2},e_{1}}$',r'$v_{q_{2},e_{2}}$'))
#plt.ylim(-0.05,1.05)


graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf",)
