#!/usr/bin/env python
'''
Script to plot various time series for Thwaites GIA experiments
'''

import sys
import os
import glob
import netCDF4
import xarray as xr
import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import time


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

# constants
rhoi = 910.0
rhow = 1028.0

def xtime2numtimeMy(xtime):
  """Define a function to convert xtime character array to numeric time values using local arithmetic"""
  # First parse the xtime character array into a string
  xtimestr = netCDF4.chartostring(xtime) # convert from the character array to an array of strings using the netCDF4 module's function

  numtime = np.zeros( (len(xtimestr),) )
  ii = 0
  dayOfMonthStart = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
  for stritem in xtimestr:
      itemarray = stritem.strip().replace('_', '-').replace(':', '-').split('-')  # Get an array of strings that are Y,M,D,h,m,s
      results = [int(i) for i in itemarray]
      numtime[ii] = results[0] + (dayOfMonthStart[results[1]-1]-1 + results[2]) / 365.0  # decimal year
      ii += 1
  return numtime

def GTtoSL(GT):
   #return (GT-VAF[0])*1.0e12/1028.0/1000.0**3/362.0e6*1000.0*1000.0
   return GT *1.0e12/1028.0/1000.0**3/362.0e6*1000.0*1000.0 * -1.0

def grounded(cellMask):
   return ((cellMask&32)//32) & ~((cellMask&4)//4)

# Define Run class for all std analyis
class modelRun:
   def __init__(self, run):
      '''
      This reads results from a model run and saves and analyzes the needed results.
      run = dictionary describing the run
      '''
      self.info = run
      self.gs = globalStats(run)
      tic = time.perf_counter()
      #self.out = outputData(run)
      toc = time.perf_counter()
      print(f"  Processed {run['name']} in {toc - tic:0.4f} seconds")

class globalStats:
   def __init__(self, run):
      # --------
      # Analysis from global stats file
      # --------
      f = netCDF4.Dataset(run['path'] + '/' + 'globalStats.nc', 'r')
      self.nt = len(f.dimensions['Time'])
      self.yrs = np.zeros((self.nt,))
      #yrs = f.variables['daysSinceStart'][:] / 365.0
      self.xtimes = f.variables['xtime'][:]


      self.yrs = xtime2numtimeMy(self.xtimes)
      #self.dyrs = self.yrs[1:] - self.yrs[0:-1]
      self.VAF = f.variables['volumeAboveFloatation'][:] / 1.0e12 * rhoi
      self.melt = f.variables['totalFloatingBasalMassBal'][:] / -1.0e12 # Gt
      self.GA = f.variables['groundedIceArea'][:] / 1000.0**2  # km^2
      self.GLflux = f.variables['groundingLineFlux'][:] / 1.0e12 # Gt/y
      self.floatArea = f.variables['floatingIceArea'][:] / 1000.0**2 # km2
      self.floatVol = f.variables['floatingIceVolume'][:] / 1000.0**3 # km3
      self.floatThk = self.floatVol / self.floatArea * 1000.0 # m
      self.grdArea = f.variables['groundedIceArea'][:] / 1000.0**2 # km2
      self.grdVol = f.variables['groundedIceVolume'][:] / 1000.0**3 # km3
      self.grdThk = self.grdVol / self.grdArea * 1000.0 # m

      self.GAloss = self.GA[0] - self.GA[:]


      # Redo on resampled time axis



class outputData:
   def __init__(self, run):
      # --------
      # Analysis from output file
      # --------
      DS = xr.open_mfdataset(run['path'] + '/' + 'output_2*.nc', combine='nested', concat_dim='Time', decode_timedelta=False, chunks={"Time": 5})
      self.yrs = DS['daysSinceStart'].load() / 365.0 + 2015.0
      self.nt = DS.dims['Time']

      cellMask = DS['cellMask']
      bed = DS['bedTopography']
      thickness = DS['thickness']
      areaCell = DS['areaCell'].load()

      # volume
      self.vol = np.zeros((self.nt,))
      self.vaf = np.zeros((self.nt,))
      for t in range(self.nt):
         self.vol[t] = (thickness.isel(Time=t) * areaCell).sum()
         self.vaf[t] = (areaCell * grounded(cellMask.isel(Time=t)) * 
            ( thickness.isel(Time=t) + (rhow / rhoi) * np.minimum(0.0*bed.isel(Time=t), bed.isel(Time=t)) )).sum()

      # uplift
      bed0 = bed.isel(time=0).load()
      self.maxUplift = np.zeros((self.nt,))
      self.maxGrdUplift = np.zeros((self.nt,))
      for t in range(self.nt):
         self.maxUplift[t] = (bed.isel(Time.t)-bed0).max()
         self.maxGrdUplift[t] = ((bed.isel(Time=t)-bed0) * grounded(cellMask.isel(Time=t))).max()


# execute analysis on each run and build output
#runData = {}  # init empty dictionary
for run in runs:
   print("Processing run: " + run['name'])

   # Build a list that contains all run data objects
   #runData[run['name']] = modelRun(run)
   run['data'] = modelRun(run)


# --------
# VAF figure
# --------
fig = plt.figure(1, facecolor='w')
axVAF = fig.add_subplot(1, 1, 1)

for run in runs:
    print("Plotting run: " + run['name'])
    yrs = run['data'].gs.yrs
    VAF = run['data'].gs.VAF

    axVAF.plot(yrs, VAF, label = run['name'], color=run['color'])

axVAF.legend(loc='best', ncol=1)
axVAF.set_xlabel('Year')
axVAF.set_ylabel('VAF (Gt)')

axSLR=axVAF.twinx()
y1, y2=axVAF.get_ylim()
x1, x2=axVAF.get_xlim()
axSLR.set_ylim(GTtoSL(y1) - GTtoSL(VAF[0]), GTtoSL(y2) - GTtoSL(VAF[0]))
#axSLR.set_yticks( range(int(GTtoSL(y1)), int(GTtoSL(y2))) )
axSLR.set_ylabel('S.L. equiv. (mm)')
axSLR.set_xlim(x1, x2)


# -------
# Bunch of stats
# -------
fig = plt.figure(2, facecolor='w', figsize=(8, 14))
nr = 6
axVAF = fig.add_subplot(nr, 1, 1)
axGA = fig.add_subplot(nr, 1, 2)
axMelt = fig.add_subplot(nr, 1, 3)
axFloatThk = fig.add_subplot(nr, 1, 4)
axGrdThk = fig.add_subplot(nr, 1, 5)
axGLf = fig.add_subplot(nr, 1, 6)

for run in runs:
    print("Plotting run: " + run['name'])
    yrs = run['data'].gs.yrs
    VAF = run['data'].gs.VAF

    axVAF.plot(yrs, run['data'].gs.VAF, label = run['name'], color=run['color'])
    axGA.plot (yrs, run['data'].gs.GA, label = run['name'], color=run['color'])
    axMelt.plot(yrs, run['data'].gs.melt, label = run['name'], color=run['color'])
    axFloatThk.plot(yrs, run['data'].gs.floatThk, label = run['name'], color=run['color'])
    axGrdThk.plot(yrs, run['data'].gs.grdThk, label = run['name'], color=run['color'])
    axGLf.plot(yrs, run['data'].gs.GLflux, label = run['name'], color=run['color'])

axVAF.legend(loc='best', ncol=1)
axVAF.set_xlabel('Year')
axVAF.set_ylabel('VAF (Gt)')

axSLR=axVAF.twinx()
y1, y2=axVAF.get_ylim()
x1, x2=axVAF.get_xlim()
axSLR.set_ylim(GTtoSL(y1) - GTtoSL(VAF[0]), GTtoSL(y2) - GTtoSL(VAF[0]))
#axSLR.set_yticks( range(int(GTtoSL(y1)), int(GTtoSL(y2))) )
axSLR.set_ylabel('S.L. equiv. (mm)')
axSLR.set_xlim(x1, x2)


#axGA.legend(loc='best', ncol=1)
axGA.set_xlabel('Year')
axGA.set_ylabel('Grounded area (km$^2$)')

axMelt.set_xlabel('Year')
axMelt.set_ylabel('Ice shelf melt rate (Gt yr${^-1}$)')
axMelt.set_ylim((0.0, 2000.0))

axFloatThk.set_xlabel('Year')
axFloatThk.set_ylabel('Mean floating ice thickness (m)')

axGrdThk.set_xlabel('Year')
axGrdThk.set_ylabel('Mean grounded ice thickness (m)')

axGLf.set_xlabel('Year')
axGLf.set_ylabel('Grounding line flux (Gt yr${^-1}$)')

# -------
# Bunch of stats vs GA
# -------
fig = plt.figure(3, facecolor='w', figsize=(8, 14))
nr = 6
axVAF = fig.add_subplot(nr, 1, 1)
axGA = fig.add_subplot(nr, 1, 2)
axMelt = fig.add_subplot(nr, 1, 3)
axFloatThk = fig.add_subplot(nr, 1, 4)
axGrdThk = fig.add_subplot(nr, 1, 5)
axGLf = fig.add_subplot(nr, 1, 6)

for run in runs:
    print("Plotting run: " + run['name'])
    GA = run['data'].gs.GAloss

    axVAF.plot(GA, run['data'].gs.VAF, label = run['name'], color=run['color'])
    axGA.plot (GA, run['data'].gs.GA, label = run['name'], color=run['color'])
    axMelt.plot(GA, run['data'].gs.melt, label = run['name'], color=run['color'])
    axFloatThk.plot(GA, run['data'].gs.floatThk, label = run['name'], color=run['color'])
    axGrdThk.plot(GA, run['data'].gs.grdThk, label = run['name'], color=run['color'])
    axGLf.plot(GA, run['data'].gs.GLflux, label = run['name'], color=run['color'])

axVAF.legend(loc='best', ncol=1)
axVAF.set_xlabel('Grounded area loss (km$^2$)')
axVAF.set_ylabel('VAF (Gt)')

axSLR=axVAF.twinx()
y1, y2=axVAF.get_ylim()
x1, x2=axVAF.get_xlim()
axSLR.set_ylim(GTtoSL(y1) - GTtoSL(VAF[0]), GTtoSL(y2) - GTtoSL(VAF[0]))
#axSLR.set_yticks( range(int(GTtoSL(y1)), int(GTtoSL(y2))) )
axSLR.set_ylabel('S.L. equiv. (mm)')
axSLR.set_xlim(x1, x2)


#axGA.legend(loc='best', ncol=1)
axGA.set_xlabel('Grounded area loss (km$^2$)')
axGA.set_ylabel('Grounded area loss (km$^2$)')

axMelt.set_xlabel('Grounded area loss (km$^2$)')
axMelt.set_ylabel('Ice shelf melt rate (Gt yr${^-1}$)')
axMelt.set_ylim((0.0, 2000.0))

axFloatThk.set_xlabel('Grounded area loss (km$^2$)')
axFloatThk.set_ylabel('Mean floating ice thickness (m)')

axGrdThk.set_xlabel('Grounded area loss (km$^2$)')
axGrdThk.set_ylabel('Mean grounded ice thickness (m)')

axGLf.set_xlabel('Grounded area loss (km$^2$)')
axGLf.set_ylabel('Grounding line flux (Gt yr${^-1}$)')


plt.show()


# --------
# uplift time series figure
# --------
fig = plt.figure('up', facecolor='w')
axUp = fig.add_subplot(1, 2, 1)
axUp2 = fig.add_subplot(1, 2, 2)
for run in runs:
    print("Plotting run: " + run['name'])
    yrs = run['data'].out.yrs
    maxUplift = run['data'].out.maxUplift
    maxGrdUplift = run['data'].out.maxGrdUplift
    vaf = run['data'].out.vaf

    axUp.plot(yrs, maxUplift, label = run['name'], color=run['color'])
    axUp.plot(yrs, maxGrdUplift, '--', color=run['color'])

    axUp2.plot(vaf[0]-vaf, maxUplift, label = run['name'], color=run['color'])
    axUp2.plot(vaf[0]-vaf, maxGrdUplift, '--', color=run['color'])

axUp.legend(loc='best', ncol=1)
axUp.set_xlabel('Year')
axUp.set_ylabel('Uplift (m)')

axUp2.legend(loc='best', ncol=1)
axUp2.set_xlabel('VAF loss (Gt)')
axUp2.set_ylabel('Uplift (m)')


# --------
plt.show()

