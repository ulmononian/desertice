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
import matplotlib.tri as tri
import cmocean

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
  "path" : "/global/cfs/cdirs/piscees/MALI_projects/Thwaites_GIA/final_analysis/control",
  "u2" : 0.0,
  "u1" : 0.0,
  "UMthk" : 0.0,
  "LT" : 0.0,
  "D" : 0.0,
  "color" : "black"
  }

N1 = {
 "name" : "N1",
 "path" : "/global/cfs/cdirs/piscees/MALI_projects/Thwaites_GIA/final_analysis/with_elastic/N1",
 "desc" : "Typical",
  "u2" : 2.0e20,
  "u1" : 4.0e21,
  "UMthk" : 670.0,
  "LT" : 60.0,
  "D" : 13.0e23,
  "color" : 'gold'
 }

N2 = {
 "name" : "N2",
 "desc" : "best2",
 "path" : "/global/cfs/cdirs/piscees/MALI_projects/Thwaites_GIA/final_analysis/with_elastic/N2",
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
 "path" : "/global/cfs/cdirs/piscees/MALI_projects/Thwaites_GIA/final_analysis/with_elastic/N3",
  "u2" : 1.0e18,
  "u1" : 1.0e18,
  "UMthk" : 0.0,
  "LT" : 25.0,
  "D" : 1.0e23,
  "color" : 'darkviolet'
 }

N4 = {
 "name" : "N4",
 "desc" : "LV-StdLith",
 "path" : "/global/cfs/cdirs/piscees/MALI_projects/Thwaites_GIA/final_analysis/with_elastic/N4",
  "u2" : 1.0e18,
  "u1" : 1.0e18,
  "UMthk" : 0.0,
  "LT" : 60.0,
  "D" : 13.0e23,
  "color" : 'firebrick'
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

# Generate triangulation to use for plotting data

def dist(i1, i2):  # helper distance fn
    #print(i1, i2, xCell[i1], xCell[i2])
    dist = ((xCell[i1]-xCell[i2])**2 + (yCell[i1]-yCell[i2])**2)**0.5
    return dist

triang = tri.Triangulation(xCell, yCell)
triMask = np.zeros(len(triang.triangles))
maxDist = 16000.0  # maximum distance in m of edges between points.  Should be somewhat larger than the coarsest mesh density (which is 8km)
for t in range(len(triang.triangles)):
    thisTri = triang.triangles[t, :]
    if dist(thisTri[0], thisTri[1]) > maxDist:
        triMask[t] = True
    if dist(thisTri[1], thisTri[2]) > maxDist:
        triMask[t] = True
    if dist(thisTri[0], thisTri[2]) > maxDist:
        triMask[t] = True
triang.set_mask(triMask)


for ax in axs.ravel(): 
    #bedPlot.append(ax.tricontourf(xCell, yCell, bed[0,:], levels=100, vmin=-2500.0, vmax=0.0, cmap='BrBG'))
    bdPlt = ax.tricontourf(triang, bed[0,:], levels=100, vmin=-2000.0, vmax=0.0, cmap=cmocean.cm.deep_r, extend='both')
    initExtentPlot.append(ax.tricontour(triang, initialExtentMask[0,:], colors='white'))
    initExtentPlot.append(ax.tricontour(triang, groundingLineMask[0,:], colors='gray'))
    ax.axis('equal')
    ax.set_xticks([])
    ax.set_xticks([], minor=True)
    ax.set_yticks([])
    ax.set_yticks([], minor=True)
f.close()

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.9, 0.1, 0.04, 0.8])
#bdPlt.set_clim(-2.0, 2.0)
cbar = fig.colorbar(bdPlt, cax=cbar_ax)


   
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
        axs[col].tricontour(triang, groundingLineMask[timeLev, :], 
           levels=[0.9999], colors=run['color'], linestyles='solid') #, label=run['name'])
    if col == 0 :
        axs[col].legend(loc='best', ncol=1)
#plt.clim(-2000,0)
        

# =============================
# Plot uplift maps

figU = plt.figure(2, facecolor='w', figsize=(8,8))
axsU = figU.subplots(2, 2)#, sharex=True, sharey=True)

run = runs[2] # best2
#f = Dataset(run['path'] + '/' + 'iceload_all.nc', 'r')
#x = f.variables['x'][:]
#y = f.variables['y'][:]
#bed = f.variables['bas']
DS = xr.open_mfdataset(run['path'] + '/' + 'output_2*.nc', combine='nested', concat_dim='Time', decode_timedelta=False, chunks={"Time": 10})
cellMask = DS["cellMask"]
#thickness = DS['thickness']
bed = DS['bedTopography']
#floatMask = (cellMask & floatValue) // floatValue
#dynamicMask = (cellMask & dynamicValue) // dynamicValue
groundingLineMask = (cellMask & groundingLineValue) // groundingLineValue
 
#axsU[0,0].pcolormesh(x, y, bed[250,:,:]-bed[0,:,:])
data=bed[250,:]-bed[0,:]
bdPlt=axsU[0,0].tricontourf(triang, data, levels=np.linspace(0.0, data.max(), num=300))
tk=np.arange(0.0, data.max(), 10.0)
#axsU[0,0].tricontour(triang, initialExtentMask[0,:]*0+1, colors='purple', levels=(0.5,))
axsU[0,0].tricontour(triang, initialExtentMask[0,:], colors='white')
axsU[0,0].tricontour(triang, groundingLineMask[0,:], colors='gray')
axsU[0,0].tricontour(triang, groundingLineMask[250,:], colors='black')
figU.colorbar(bdPlt, ax=axsU[0,0], ticks=tk)
axsU[0,0].set_xticks([])
axsU[0,0].set_xticks([], minor=True)
axsU[0,0].set_yticks([])
axsU[0,0].set_yticks([], minor=True)
axsU[0,0].annotate('a', (0.03, 0.95), xycoords='axes fraction', fontweight='bold', fontsize='large')

#axsU[0,1].pcolormesh(x, y, bed[500,:,:]-bed[0,:,:])
#axsU[0,0].tricontour(triang, initialExtentMask[0,:]*0+1, colors='purple', levels=(0.5,))
data=bed[500,:]-bed[0,:]
bdPlt=axsU[0,1].tricontourf(triang, data, levels=np.linspace(0.0, data.max(), num=300))
axsU[0,1].tricontour(triang, initialExtentMask[0,:], colors='white')
axsU[0,1].tricontour(triang, groundingLineMask[0,:], colors='gray')
axsU[0,1].tricontour(triang, groundingLineMask[500,:], colors='black')
tk=np.arange(0.0, data.max(), 20.0)
figU.colorbar(bdPlt, ax=axsU[0,1], ticks=tk)
axsU[0,1].set_xticks([])
axsU[0,1].set_xticks([], minor=True)
axsU[0,1].set_yticks([])
axsU[0,1].set_yticks([], minor=True)
axsU[0,1].annotate('b', (0.03, 0.95), xycoords='axes fraction', fontweight='bold', fontsize='large')
     
     
run = runs[3] # best2
DS = xr.open_mfdataset(run['path'] + '/' + 'output_2*.nc', combine='nested', concat_dim='Time', decode_timedelta=False, chunks={"Time": 10})
bed = DS['bedTopography']
cellMask = DS["cellMask"]
groundingLineMask = (cellMask & groundingLineValue) // groundingLineValue
#f = Dataset(run['path'] + '/' + 'iceload_all.nc', 'r')
#x = f.variables['x'][:]
#y = f.variables['y'][:]
#bed = f.variables['bas']

#axsU[0,0].pcolormesh(x, y, bed[250,:,:]-bed[0,:,:])
data=bed[250,:]-bed[0,:]
bdPlt=axsU[1,0].tricontourf(triang, data, levels=np.linspace(0.0, data.max(), num=300))
#axsU[0,0].tricontour(triang, initialExtentMask[0,:]*0+1, colors='purple', levels=(0.5,))
axsU[1,0].tricontour(triang, initialExtentMask[0,:], colors='white')
axsU[1,0].tricontour(triang, groundingLineMask[0,:], colors='gray')
axsU[1,0].tricontour(triang, groundingLineMask[250,:], colors='black')
tk=np.arange(0.0, data.max(), 10.0)
figU.colorbar(bdPlt, ax=axsU[1,0], ticks=tk)
axsU[1,0].set_xticks([])
axsU[1,0].set_xticks([], minor=True)
axsU[1,0].set_yticks([])
axsU[1,0].set_yticks([], minor=True)
axsU[1,0].annotate('c', (0.03, 0.95), xycoords='axes fraction', fontweight='bold', fontsize='large')

#axsU[0,1].pcolormesh(x, y, bed[500,:,:]-bed[0,:,:])
#axsU[0,0].tricontour(triang, initialExtentMask[0,:]*0+1, colors='purple', levels=(0.5,))
data=bed[500,:]-bed[0,:]
bdPlt=axsU[1,1].tricontourf(triang, data, levels=np.linspace(0.0, data.max(), num=300))
axsU[1,1].tricontour(triang, initialExtentMask[0,:], colors='white')
axsU[1,1].tricontour(triang, groundingLineMask[0,:], colors='gray')
axsU[1,1].tricontour(triang, groundingLineMask[500,:], colors='black')
tk=np.arange(0.0, data.max(), 20.0)
figU.colorbar(bdPlt, ax=axsU[1,1], ticks=tk)
axsU[1,1].set_xticks([])
axsU[1,1].set_xticks([], minor=True)
axsU[1,1].set_yticks([])
axsU[1,1].set_yticks([], minor=True)
axsU[1,1].annotate('d', (0.03, 0.95), xycoords='axes fraction', fontweight='bold', fontsize='large')
    
# Customize plots
#for ax in axs:
#   ax.axis('equal')
#axs[0].legend(loc='best', ncol=1)
#axs[0,0].set_ylim(top=-1020000, bottom = -1120000) # hard-code limits for now
#axs[0,0].set_xlim(left=-425000, right=-300000)



plt.show()
