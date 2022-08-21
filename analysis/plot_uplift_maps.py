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
from matplotlib.colors import TwoSlopeNorm

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
timeLevs = [85, 285, 500]
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

f = Dataset(runs[0]['path']+"/output_2015.nc", 'r')
xCell = f.variables["xCell"][:]
yCell = f.variables["yCell"][:]
cellMask = f.variables["cellMask"][:]
thickness0 = f.variables["thickness"][0]
initialExtentMask = (cellMask & initialExtentValue) // initialExtentValue
initialExtentMask2 = (thickness0>2000.0)
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

boundary = set()
for i in range(len(triang.neighbors)):
    for j in range(3):
        if (triang.neighbors[i,j] == -1):
            boundary.add(triang.triangles[i,j])
            boundary.add(triang.triangles[i,(j+1)%3])
#            nk1,nk2 = (k+1)%3, (k+2)%3
#            boundary.add(triang.simplices[i][nk1])
#            boundary.add(triang.simplices[i][nk2])

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

# GIA file
#DSgia = xr.open_mfdataset(run['path'] + '/' + 'uplift_GIA_all.nc', combine='nested', concat_dim='Time', decode_timedelta=False, chunks={"Time": 10})
DSgia = Dataset(run['path'] + '/' + 'uplift_GIA_all.nc')
uplift = DSgia.variables['uplift']
xgia = DSgia.variables['x'][:]
ygia = DSgia.variables['y'][:]

tk2=np.arange(0.0, 268.0, 20.0)
levs=np.arange(0.0, 268.0, 0.1)

#axsU[0,0].pcolormesh(x, y, bed[250,:,:]-bed[0,:,:])
data=bed[285,:]-bed[0,:]
#bdPlt=axsU[0,0].tricontourf(triang, data, levels=np.linspace(0.0, data.max(), num=300))
data = uplift[285][:,:]
norm = TwoSlopeNorm(vmin=np.amin(data), vcenter=0, vmax=np.amax(data))
#bdPlt=axsU[0,0].contourf(xgia, ygia, data, levels=3000, norm=norm, cmap='RdBu') # include subsidence
bdPlt=axsU[0,0].contourf(xgia, ygia, data, levels=levs, extend='min', cmap='viridis_r')
axsU[0,0].contour(xgia, ygia, data, levels=np.arange(-10, 0.2, 1.0), colors='olive')
tk=np.arange(0.0, data.max(), 10.0)
#tk=np.insert(tk[1:],0,np.amin(data))
#axsU[0,0].tricontour(triang, initialExtentMask[0,:]*0+1, colors='purple', levels=(0.5,))
axsU[0,0].tricontour(triang, initialExtentMask[0,:], colors='white')
#axsU[0,0].tricontour(triang, initialExtentMask2, colors='white', levels=(0.5,))
axsU[0,0].tricontour(triang, groundingLineMask[0,:], colors='gray')
axsU[0,0].tricontour(triang, groundingLineMask[285,:], colors='black')
axsU[0,0].plot(xCell[list(boundary)], yCell[list(boundary)], '.r', markersize=1) # plot MALI domain bdy
#print(f'shape={data.shape}')
#print(f'len xgia={xgia.shape}, len ygia={ygia.shape}')
(yind,xind)=np.unravel_index(np.argmax(data), data.shape)
axsU[0,0].plot(xgia[xind], ygia[yind], 'wx')
axsU[0,0].annotate(f"{data[yind,xind]:.2f} m", (xgia[xind], ygia[yind]), color='w', textcoords='offset points', xytext=(5,5))
#print(f'max={data.max()}, xyval={data[yind,xind]}, xyval2={data[xind,yind]}')
cbar = figU.colorbar(bdPlt, ax=axsU[0,0], ticks=tk2)
cbar.set_label('uplift (m)', rotation=270, labelpad=15)
axsU[0,0].set_xticks([])
axsU[0,0].set_xticks([], minor=True)
axsU[0,0].set_yticks([])
axsU[0,0].set_yticks([], minor=True)
axsU[0,0].annotate('a', (0.03, 0.95), xycoords='axes fraction', fontweight='bold', fontsize='large')

#axsU[0,1].pcolormesh(x, y, bed[500,:,:]-bed[0,:,:])
#axsU[0,0].tricontour(triang, initialExtentMask[0,:]*0+1, colors='purple', levels=(0.5,))
data=bed[500,:]-bed[0,:]
#bdPlt=axsU[0,1].tricontourf(triang, data, levels=np.linspace(0.0, data.max(), num=300))
data = uplift[500][:,:]
norm = TwoSlopeNorm(vmin=np.amin(data), vcenter=0, vmax=np.amax(data))
#bdPlt=axsU[0,1].contourf(xgia, ygia, data, levels=3000, norm=norm, cmap='RdBu')
bdPlt=axsU[0,1].contourf(xgia, ygia, data, levels=levs, extend='min', cmap='viridis_r')
axsU[0,1].contour(xgia, ygia, data, levels=np.arange(-10, 0.2, 1.0), colors='olive')
axsU[0,1].tricontour(triang, initialExtentMask[0,:], colors='white')
axsU[0,1].tricontour(triang, groundingLineMask[0,:], colors='gray')
axsU[0,1].tricontour(triang, groundingLineMask[500,:], colors='black')
tk=np.arange(0.0, data.max(), 20.0)
#tk=np.insert(tk[1:],0,np.amin(data))
axsU[0,1].plot(xCell[list(boundary)], yCell[list(boundary)], '.r', markersize=1) # plot MALI domain bdy
(yind,xind)=np.unravel_index(np.argmax(data), data.shape)
axsU[0,1].plot(xgia[xind], ygia[yind], 'wx')
axsU[0,1].annotate(f"{data[yind,xind]:.2f} m", (xgia[xind], ygia[yind]), color='w', textcoords='offset points', xytext=(5,5))
cbar = figU.colorbar(bdPlt, ax=axsU[0,1], ticks=tk2)
cbar.set_label('uplift (m)', rotation=270, labelpad=15)
axsU[0,1].set_xticks([])
axsU[0,1].set_xticks([], minor=True)
axsU[0,1].set_yticks([])
axsU[0,1].set_yticks([], minor=True)
axsU[0,1].annotate('b', (0.03, 0.95), xycoords='axes fraction', fontweight='bold', fontsize='large')
     
     
run = runs[3] # VLV-thin
DS = xr.open_mfdataset(run['path'] + '/' + 'output_2*.nc', combine='nested', concat_dim='Time', decode_timedelta=False, chunks={"Time": 10})
bed = DS['bedTopography']
cellMask = DS["cellMask"]
groundingLineMask = (cellMask & groundingLineValue) // groundingLineValue
#f = Dataset(run['path'] + '/' + 'iceload_all.nc', 'r')
#x = f.variables['x'][:]
#y = f.variables['y'][:]
#bed = f.variables['bas']

# GIA file
#DSgia = xr.open_mfdataset(run['path'] + '/' + 'uplift_GIA_all.nc', combine='nested', concat_dim='Time', decode_timedelta=False, chunks={"Time": 10})
#uplift = DSgia['uplift']
#xgia = DSgia['x']
#ygia = DSgia['y']
DSgia = Dataset(run['path'] + '/' + 'uplift_GIA_all.nc')
uplift = DSgia.variables['uplift']
xgia = DSgia.variables['x'][:]
ygia = DSgia.variables['y'][:]

#axsU[0,0].pcolormesh(x, y, bed[250,:,:]-bed[0,:,:])
data=bed[285,:]-bed[0,:]
#bdPlt=axsU[1,0].tricontourf(triang, data, levels=np.linspace(0.0, data.max(), num=300))
data = uplift[285][:,:]
norm = TwoSlopeNorm(vmin=np.amin(data), vcenter=0, vmax=np.amax(data))
#bdPlt=axsU[1,0].contourf(xgia, ygia, data, levels=3000, norm=norm, cmap='RdBu')
bdPlt=axsU[1,0].contourf(xgia, ygia, data, levels=levs, extend='min', cmap='viridis_r')
axsU[1,0].contour(xgia, ygia, data, levels=np.arange(-10, 0.2, 1.0), colors='olive')
#axsU[0,0].tricontour(triang, initialExtentMask[0,:]*0+1, colors='purple', levels=(0.5,))
axsU[1,0].tricontour(triang, initialExtentMask[0,:], colors='white')
axsU[1,0].tricontour(triang, groundingLineMask[0,:], colors='gray')
axsU[1,0].tricontour(triang, groundingLineMask[285,:], colors='black')
tk=np.arange(0.0, data.max(), 10.0)
#tk=np.insert(tk[1:],0,np.amin(data))
axsU[1,0].plot(xCell[list(boundary)], yCell[list(boundary)], '.r', markersize=1) # plot MALI domain bdy
(yind,xind)=np.unravel_index(np.argmax(data), data.shape)
axsU[1,0].plot(xgia[xind], ygia[yind], 'wx')
axsU[1,0].annotate(f"{data[yind,xind]:.2f} m", (xgia[xind], ygia[yind]), color='w', textcoords='offset points', xytext=(5,5))
cbar = figU.colorbar(bdPlt, ax=axsU[1,0], ticks=tk2)
cbar.set_label('uplift (m)', rotation=270, labelpad=15)
axsU[1,0].set_xticks([])
axsU[1,0].set_xticks([], minor=True)
axsU[1,0].set_yticks([])
axsU[1,0].set_yticks([], minor=True)
axsU[1,0].annotate('c', (0.03, 0.95), xycoords='axes fraction', fontweight='bold', fontsize='large')

#axsU[0,1].pcolormesh(x, y, bed[500,:,:]-bed[0,:,:])
#axsU[0,0].tricontour(triang, initialExtentMask[0,:]*0+1, colors='purple', levels=(0.5,))
data=bed[500,:]-bed[0,:]
#bdPlt=axsU[1,1].tricontourf(triang, data, levels=np.linspace(0.0, data.max(), num=300))
data = uplift[500][:,:]
norm = TwoSlopeNorm(vmin=np.amin(data), vcenter=0, vmax=np.amax(data))
#bdPlt=axsU[1,1].contourf(xgia, ygia, data, levels=3000, norm=norm, cmap='RdBu')
bdPlt=axsU[1,1].contourf(xgia, ygia, data, levels=levs, extend='min', cmap='viridis_r')
axsU[1,1].contour(xgia, ygia, data, levels=np.arange(-10, 0.2, 1.0), colors='olive')
#bdPlt=axsU[1,1].contourf(xgia, ygia, uplift[500], levels=300)
axsU[1,1].tricontour(triang, initialExtentMask[0,:], colors='white')
axsU[1,1].tricontour(triang, groundingLineMask[0,:], colors='gray')
axsU[1,1].tricontour(triang, groundingLineMask[500,:], colors='black')
tk=np.arange(0.0, data.max(), 20.0)
#tk=np.insert(tk[1:],0,np.amin(data))
axsU[1,1].plot(xCell[list(boundary)], yCell[list(boundary)], '.r', markersize=1) # plot MALI domain bdy
(yind,xind)=np.unravel_index(np.argmax(data), data.shape)
axsU[1,1].plot(xgia[xind], ygia[yind], 'wx')
axsU[1,1].annotate(f"{data[yind,xind]:.2f} m", (xgia[xind], ygia[yind]), color='w', textcoords='offset points', xytext=(5,5))
cbar = figU.colorbar(bdPlt, ax=axsU[1,1], ticks=tk2)
cbar.set_label('uplift (m)', rotation=270, labelpad=15)
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
