#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys, math, decimal, os.path
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
from matplotlib.ticker import AutoMinorLocator
from matplotlib import rcParams
plt.style.use('base')


rcParams['text.usetex'] = True 
rcParams['font.family'] = 'sans-serif'
rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]

file = sys.argv[1]

f = open(file)
fl = f.readlines()
f.close()

len_file = len(fl)

data = pd.read_csv(file,sep=";",nrows=500000)

print(data.describe())
if "time" not in data.columns.values:
    data["time"] = range(0,data.shape[0])


subset = data[(data["time"] > (50000 - 10)) & (data["time"] < data["time"].max() - 10)]


subset["time"] = subset["time"] - (50000 - 10)


def log_the_time(row):
    return(math.log10(float(row["time"])))

subset["logtime"] = subset.apply(log_the_time,axis=1)

fig = plt.figure(figsize=(7.5,7.5))
gs = gridspec.GridSpec(1,1)

ax = plt.subplot(gs[0,0])
ax.plot(subset["logtime"], subset["mbar"], color="red", linewidth=2)

graphname = "graph_" + os.path.splitext(os.path.basename(file))[0]
plt.savefig(graphname,format="pdf", bbox_inches="tight")
