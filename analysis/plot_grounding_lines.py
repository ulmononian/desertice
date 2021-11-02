#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 15:40:48 2021

@author: trevorhillebrand
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import os
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import xarray as xr
from optparse import OptionParser
import matplotlib.pyplot as plt
import matplotlib as mpl

rhoi = 910.0
rhosw = 1028.0
#print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
#parser = OptionParser(description=__doc__)
#parser.add_option("-r", dest="runs", help="path to .nc file or dir containing output.nc file (strings separated by commas; no spaces)", default=None, metavar="FILENAME")
#parser.add_option("-t", dest="timeLevels", help="integer time levels at which to plot (int separated by commas; no spaces)", default=-1, metavar="FILENAME")
#
#options, args = parser.parse_args()
##if type(options.runs) is list: 
#runs = options.runs.split(',') # split run directories into list
##else:
##    runs = [options.runs] # Kind of cludgey, but allows for consistent indexing in loop
    
#if type(options.timeLevels) is list:
#timeLevs = options.timeLevels.split(',') # split time levels into list
#else:
timeLevs = [100, 250, 500]
# Define runs

ctrl = {
  "name" : "ctrl",
 "desc" : "control",
  "path" : "/global/project/projectdirs/piscees/MALI_projects/Thwaites_GIA/Thwaites_1km_GIA_ensemble_March2021/control_medianAIS",
  "u2" : 0.0,
  "u1" : 0.0,
  "UMthk" : 0.0,
  "LT" : 0.0,
  "D" : 0.0,
  "color" : "black"
  }

N1 = {
 "name" : "N1",
 "desc" : "Typical",
# "path" : "/global/cscratch1/sd/trhille/Thwaites_1km_GIA_ensemble/N1_01yr/",
 "path" : "/global/project/projectdirs/piscees/MALI_projects/Thwaites_GIA/Thwaites_1km_GIA_ensemble_March2021/N1",
  "u2" : 2.0e20,
  "u1" : 4.0e21,
  "UMthk" : 670.0,
  "LT" : 60.0,
  "D" : 13.0e23,
  "color" : 'firebrick'
 }

N2 = {
 "name" : "N2",
 "desc" : "best2",
# "path" : "/global/cscratch1/sd/trhille/Thwaites_1km_GIA_ensemble/N2_01yr/",
 "path" : "/global/project/projectdirs/piscees/MALI_projects/Thwaites_GIA/Thwaites_1km_GIA_ensemble_March2021/N2",
  "u2" : 4.0e18,
  "u1" : 2.0e19,
  "UMthk" : 200.0,
  "LT" : 60.0,
  "D" : 13.0e23,
  "color" : 'tab:orange' 
 }

N3 = {
 "name" : "N3",
 "desc" : "LV-ThinLith",
# "path" : "/global/cscratch1/sd/trhille/Thwaites_1km_GIA_ensemble/N3_1yr/",
 "path" : "/global/project/projectdirs/piscees/MALI_projects/Thwaites_GIA/Thwaites_1km_GIA_ensemble_March2021/N3",
  "u2" : 1.0e18,
  "u1" : 1.0e18,
  "UMthk" : 0.0,
  "LT" : 25.0,
  "D" : 1.0e23,
  "color" : 'gold'
 }

N4 = {
 "name" : "N4",
 "desc" : "LV-StdLith",
# "path" : "/global/cscratch1/sd/trhille/Thwaites_1km_GIA_ensemble/N4_1yr/",
 "path" : "/global/project/projectdirs/piscees/MALI_projects/Thwaites_GIA/Thwaites_1km_GIA_ensemble_March2021/N4",
  "u2" : 1.0e18,
  "u1" : 1.0e18,
  "UMthk" : 0.0,
  "LT" : 60.0,
  "D" : 13.0e23,
  "color" : 'yellowgreen'
 }

runs = [ctrl, N1, N2, N3, N4]


# Define cellMask bit values so we don't have to use bitmask conversion script
initialExtentValue = 1
dynamicValue = 2
floatValue = 4
groundingLineValue = 256

# Use rows for q values, columns for climate forcings.
# Currently this means plotting all time levels on the same
# axes, which may or may not work, depending on application.
# This could definitely be improved.
fig = plt.figure(1, facecolor='w', figsize=(16,7))
axs = fig.subplots(1, 3)#, sharex=True, sharey=True) 


# Plot bed topo and initial extent on all axes using first file
bedPlot = []
initExtentPlot = []
f = Dataset(runs[0]['path']+"/output_2015.nc", 'r')
xCell = f.variables["xCell"][:]
yCell = f.variables["yCell"][:]
bed = f.variables["bedTopography"][:]
cellMask = f.variables["cellMask"][:]
initialExtentMask = (cellMask & initialExtentValue) // initialExtentValue
groundingLineMask = (cellMask & groundingLineValue) // groundingLineValue

for ax in axs.ravel(): 
    bedPlot.append(ax.tricontourf(xCell, yCell, bed[0,:], levels=100, vmin=-2500.0, vmax=0.0, cmap='BrBG'))
    initExtentPlot.append(ax.tricontour(xCell, yCell, initialExtentMask[0,:], colors='black'))
    initExtentPlot.append(ax.tricontour(xCell, yCell, groundingLineMask[0,:], colors='gray'))
    ax.axis('equal')
f.close()

# Functions to simplify plotting
def get_line_color(run):
    if 'm5_' in run or 'm5/' in run:
        lowCalving = 'VM180'
        medCalving = 'VM170'
        highCalving = 'VM160'
    elif 'm7_' in run or 'm7/' in run:
        if 'HadGEM2' in run or 'CNRM' in run:
            lowCalving = 'VM190'
            medCalving = 'VM180'
            highCalving = 'VM170'
        elif 'MIROC5' in run:
            lowCalving = 'VM180'
            medCalving = 'VM170'
            highCalving = 'VM160'
    
    if '2017calvingFront' in run or 'calvingVelocityData' in run:
        lineColor = 'tab:grey'
    elif highCalving in run:
        lineColor = 'tab:purple'
    elif medCalving in run:
        lineColor = 'tab:blue'
    elif lowCalving in run:
        lineColor = 'tab:cyan'
    else:
        lineColor = 'white'
    
    return lineColor

def get_line_style(run):
    if 'm5_' in run or 'm5/' in run:
        lineStyle = 'solid'
    elif 'm7_' in run or 'm7/' in run:
        lineStyle = 'dashed'
    else:
        lineStyle = 'none'
        
    return lineStyle
    
#Loop over runs
for run in runs:
    DS = xr.open_mfdataset(run['path'] + '/' + 'output_2*.nc', combine='nested', concat_dim='Time', decode_timedelta=False, chunks={"Time": 10})
    cellMask = DS["cellMask"]
    thickness = DS['thickness']
    bed = DS['bedTopography']
    floatMask = (cellMask & floatValue) // floatValue
    dynamicMask = (cellMask & dynamicValue) // dynamicValue
    groundingLineMask = (cellMask & groundingLineValue) // groundingLineValue
    col = 0
    row = 0
    
    # Loop over time levels
    for t in range(len(timeLevs)):
        timeLev = int(timeLevs[t]) #these are strings for some reason; make int to index
        col = t
        axs[col].tricontour(xCell, yCell, groundingLineMask[timeLev, :], 
           levels=[0.9999], colors=run['color'], linestyles='solid') #, label=run['name'])
        

        
# Customize plots
#for ax in axs:
#   ax.axis('equal')
#axs[0].legend(loc='best', ncol=1)
#axs[0,0].set_ylim(top=-1020000, bottom = -1120000) # hard-code limits for now
#axs[0,0].set_xlim(left=-425000, right=-300000)

plt.show()
