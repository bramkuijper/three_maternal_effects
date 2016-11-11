#!/usr/bin/env python3

# plot the battleground between offspring and parent strategies

import pandas as pd
import itertools
import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import rcParams
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from matplotlib import cm

plt.style.use('base')

rcParams['axes.labelsize'] = 15
rcParams['text.usetex'] = True
rcParams['font.family'] = 'sans-serif'

# see http://stackoverflow.com/questions/2537868/sans-serif-math-with-latex-in-matplotlib 
rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]  

# initialize the figure
fig = plt.figure(figsize=(14,5))

widths = [ 1, 1, 1, 0.05 ]
heights = [ 1, 0.1]

gs = gridspec.GridSpec(len(heights), len(widths), width_ratios=widths, height_ratios=heights)

#data = pd.read_csv("../battleground_k_model/summary_battleground_k_model.csv", sep=";")
data = pd.read_csv("../resolution_k_model/summary_nonlocal_n1.csv", sep=";")

c1val = 0.83

# get a subset to do some line drawing
subset = data[
                (data["n"] == 1.001) 
                & (data["l"] == 0.5) 
                & (data["gam"] == 1.0) 
                & (data["d"] == 0.1) 
                & (data["sigma21"] == 0.25) 
                & (data["sigma12"] == 0.25)
                & (data["c1"] == c1val)
                ]

# first subplot: battleground envt e1 when gamma = 1
ax = plt.subplot(gs[0,0])

p1 = ax.plot(
                subset["c2"],
                subset["z1_e1_off"],
                color="blue",
                linewidth=2,
                label="Offspring, $e_{1}$")
p2 = ax.plot(
                subset["c2"],
                subset["z1_e1_mom"],
                color="blue",
                linestyle=(0,(3,2)),
                linewidth=2,
                label="Mother, $e_{1}$")

p3 = ax.plot(
                subset["c2"],
                subset["z1_e2_off"],
                color="red",
                linewidth=1,
                label="Offspring, $e_{2}$")
p4 = ax.plot(
                subset["c2"],
                subset["z1_e2_mom"],
                color="red",
                linestyle=(0,(3,2)),
                linewidth=1,
                label="Mother, $e_{2}$")

ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.set_ylim((-0.05,1.05))
ax.set_xlim((-0.05,1.05))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position("left")
ax.xaxis.set_ticks_position("bottom")

ax.set_ylabel(ylabel="Proportion $z_{1}$ offspring in envt $e_{i}$, $f_{i}$")

ax.set_title(label="Offspring favor more extreme mixtures", position=(0.5,1.05))
ax.set_title(loc="left", label="A", position=(0.05,0.05))

data = pd.read_csv("../resolution_k_model/summary_gamma_beta.csv", sep=";")

# get a subset to do some line drawing
subset = data[
                (data["n"] == 1.001) 
                & (data["l"] == 0.5) 
                & (data["gam1"] == 2.0) 
                & (data["beta2"] == 1.0) 
                & (data["d"] == 0.1) 
                & (data["sigma21"] == 0.25) 
                & (data["sigma12"] == 0.25)
                & (data["c1"] == c1val)
                ]

ax = plt.subplot(gs[0,1])

p1 = ax.plot(
                subset["c2"],
                subset["z1_e1_off"],
                color="blue",
                linewidth=2,
                label="Offspring, $e_{1}$")
p2 = ax.plot(
                subset["c2"],
                subset["z1_e1_mom"],
                color="blue",
                linestyle=(0,(3,2)),
                linewidth=2,
                label="Mother, $e_{1}$")

p3 = ax.plot(
                subset["c2"],
                subset["z1_e2_off"],
                color="red",
                linewidth=1,
                label="Offspring, $e_{2}$")
p4 = ax.plot(
                subset["c2"],
                subset["z1_e2_mom"],
                color="red",
                linestyle=(0,(3,2)),
                linewidth=1,
                label="Mother, $e_{2}$")


ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.set_ylim((-0.05,1.05))
ax.set_xlim((-0.05,1.05))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position("left")
ax.xaxis.set_ticks_position("bottom")
ax.set_yticklabels([])
ax.set_xlabel(xlabel=
    r"Cost of maladaptation in environment $e_{\thinspace 2}$, $c_{\thinspace 2}$")
ax.set_title(label="Offspring favor more $z_{2}$", position=(0.5,1.04))
ax.set_title(loc="left", label="B", position=(0.05,0.05))




data = pd.read_csv("../resolution_k_model/summary_gamma_beta.csv", sep=";")

# get a subset to do some line drawing
subset = data[
                (data["n"] == 1.001) 
                & (data["l"] == 0.5) 
                & (data["gam1"] == 2.0) 
                & (data["beta2"] == 2.0) 
                & (data["d"] == 0.1) 
                & (data["sigma21"] == 0.25) 
                & (data["sigma12"] == 0.25)
                & (data["c1"] == c1val)
                ]

ax = plt.subplot(gs[0,2])

p1 = ax.plot(
                subset["c2"],
                subset["z1_e1_off"],
                color="blue",
                linewidth=2,
                label="Offspring in $e_{1}$")
p2 = ax.plot(
                subset["c2"],
                subset["z1_e1_mom"],
                color="blue",
                linestyle=(0,(3,2)),
                linewidth=2,
                label="Mother in $e_{1}$")

p3 = ax.plot(
                subset["c2"],
                subset["z1_e2_off"],
                color="red",
                linewidth=1,
                label="Offspring in $e_{2}$")
p4 = ax.plot(
                subset["c2"],
                subset["z1_e2_mom"],
                color="red",
                linestyle=(0,(3,2)),
                linewidth=1,
                label="Mother in $e_{2}$")


ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.set_ylim((-0.05,1.05))
ax.set_xlim((-0.05,1.05))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position("left")
ax.xaxis.set_ticks_position("bottom")
ax.set_yticklabels([])
ax.set_title(label="Offspring favor less extreme mixtures", position=(0.5,1.05))
ax.set_title(loc="left", label="C", position=(0.05,0.05))


plt.legend(bbox_to_anchor=(1.55,1))
format = "pdf"
plt.savefig("battleground_lineplot." + format, format=format)

