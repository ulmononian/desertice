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

doGIAfiles = True
doGIAfiles = False

doUpliftTS = False


# Define runs

ctrl = {
  "name" : "ctrl",
 "desc" : "control",
  "legname" : "CTRL",
  "path" : "/global/cfs/cdirs/piscees/MALI_projects/Thwaites_GIA/final_analysis/control",
  "u2" : 0.0,
  "u1" : 0.0,
  "UMthk" : 0.0,
  "LT" : 0.0,
  "D" : 0.0,
  "color" : "black",
  "style" : "-"
  }

N1 = {
 "name" : "N1",
 "desc" : "Typical",
  "legname" : "TYP",
# "path" : "/global/cscratch1/sd/trhille/Thwaites_1km_GIA_ensemble/N1_01yr/",
 "path" : "/global/cfs/cdirs/piscees/MALI_projects/Thwaites_GIA/final_analysis/with_elastic/N1",
  "u2" : 2.0e20,
  "u1" : 4.0e21,
  "UMthk" : 670.0,
  "LT" : 60.0,
  "D" : 13.0e23,
  #"color" : 'firebrick'
  "color" : 'gold',
  "style" : "-"
 }

N2 = {
 "name" : "N2",
 "desc" : "best2",
  "legname" : "BEST2",
# "path" : "/global/cscratch1/sd/trhille/Thwaites_1km_GIA_ensemble/N2_01yr/",
 "path" : "/global/cfs/cdirs/piscees/MALI_projects/Thwaites_GIA/final_analysis/with_elastic/N2",
  "u2" : 4.0e18,
  "u1" : 2.0e19,
  "UMthk" : 200.0,
  "LT" : 60.0,
  "D" : 13.0e23,
  #"color" : 'tab:orange' 
  "color" : 'tab:orange',
  "style" : "-"
 }

N3 = {
 "name" : "N3",
 "desc" : "LV-ThinLith",
  "legname" : "VLV-THIN",
# "path" : "/global/cscratch1/sd/trhille/Thwaites_1km_GIA_ensemble/N3_1yr/",
 "path" : "/global/cfs/cdirs/piscees/MALI_projects/Thwaites_GIA/final_analysis/with_elastic/N3",
  "u2" : 1.0e18,
  "u1" : 1.0e18,
  "UMthk" : 0.0,
  "LT" : 25.0,
  "D" : 1.0e23,
  #"color" : 'gold',
  "color" : 'darkviolet',
  "style" : "-",
  "bad_uplift_ind" : np.arange(253, 265+1)
 }

N4 = {
 "name" : "N4",
 "desc" : "LV-StdLith",
  "legname" : "VLV",
# "path" : "/global/cscratch1/sd/trhille/Thwaites_1km_GIA_ensemble/N4_1yr/",
 "path" : "/global/cfs/cdirs/piscees/MALI_projects/Thwaites_GIA/final_analysis/with_elastic/N4",
  "u2" : 1.0e18,
  "u1" : 1.0e18,
  "UMthk" : 0.0,
  "LT" : 60.0,
  "D" : 13.0e23,
  #"color" : 'yellowgreen',
  "color" : 'firebrick',
  "style" : "-",
  "bad_uplift_ind" : np.array([233, 234, 235])
 }

PIGLctrl = {
  "name" : "PIGLctrl",
 "desc" : "PIGLcontrol",
  "legname" : "HM-CTRL",
  "path" : "/global/cfs/cdirs/piscees/MALI_projects/Thwaites_GIA/final_analysis/PIGL_control",
  "u2" : 0.0,
  "u1" : 0.0,
  "UMthk" : 0.0,
  "LT" : 0.0,
  "D" : 0.0,
  "color" : "gray",
  "style" : "-"
  }

PIGLN3 = {
 "name" : "PIGLN3",
 "desc" : "PIGL-LV-ThinLith",
  "legname" : "HM-VLV-THIN",
# "path" : "/global/cscratch1/sd/trhille/Thwaites_1km_GIA_ensemble/N3_1yr/",
 "path" : "/global/cfs/cdirs/piscees/MALI_projects/Thwaites_GIA/final_analysis/with_elastic/PIGL_N3",
  "u2" : 1.0e18,
  "u1" : 1.0e18,
  "UMthk" : 0.0,
  "LT" : 25.0,
  "D" : 1.0e23,
  #"color" : 'deepskyblue'
  "color" : 'plum',
  "style" : "-"
 }


runs = [ctrl, N1, N2, N4, N3]  # standard set of runs
##runs = [ctrl, N1, N2] # testing faster
#runs = [ctrl, PIGLctrl, N3, PIGLN3]; doGIAfiles = False  # melt comparison runs

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
      tic = time.perf_counter()
      self.info = run
      self.gs = globalStats(run)
      if doGIAfiles:
         self.GIA = GIAoutputData(run)
      if doUpliftTS:
         self.out = outputData(run)
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
      # Clean a few melt blips
      self.melt[self.melt>3000.0] = np.NAN
      self.meltrate = f.variables['avgSubshelfMelt'][:]
      self.meltrate[self.meltrate>500.0] = np.NAN
      self.GA = f.variables['groundedIceArea'][:] / 1000.0**2  # km^2
      self.GLflux = f.variables['groundingLineFlux'][:] / 1.0e12 # Gt/y
      self.GLflux[0] = np.NAN # remove the zero in the first time level
      self.floatArea = f.variables['floatingIceArea'][:] / 1000.0**2 # km2
      self.floatVol = f.variables['floatingIceVolume'][:] / 1000.0**3 # km3
      self.floatThk = self.floatVol / self.floatArea * 1000.0 # m
      self.grdArea = f.variables['groundedIceArea'][:] / 1000.0**2 # km2
      self.grdVol = f.variables['groundedIceVolume'][:] / 1000.0**3 # km3
      self.grdThk = self.grdVol / self.grdArea * 1000.0 # m

      self.VAFrate = np.zeros(self.VAF.shape)
      self.VAFrate[1:] = (self.VAF[1:] - self.VAF[:-1]) / (self.yrs[1:] - self.yrs[:-1])

      self.GAloss = self.GA[0] - self.GA[:]

      # Redo on resampled time axis
      resampEndtime = 500.0
      self.resampYrs = np.linspace(2015.0, 2015.0+resampEndtime, num=int(resampEndtime*12*1))
      self.resampVAF = np.interp(self.resampYrs, self.yrs, self.VAF) # generate y values for each x
      self.resampVAFrate = (self.resampVAF[1:] - self.resampVAF[0:-1]) / (self.resampYrs[1:] - self.resampYrs[0:-1])
      self.resampVAFrate[self.resampVAFrate<6*self.resampVAFrate.mean()] = np.nan # remove some outliers
      # remove more outliers
      for i in range(len(self.resampVAFrate)):
          if self.resampVAFrate[i] < 2.0 * self.resampVAFrate[i-20:i+20].mean():
              self.resampVAFrate[i] = np.nan


      # calculate SLR reduction
      if not run['name'] in  ('ctrl', 'PIGLctrl'):
          ctrlVAF = runs[0]['data'].gs.resampVAF
          self.reduction = (1.0 - (self.resampVAF - self.resampVAF[0]) / (ctrlVAF - ctrlVAF[0])) * 100.0

      # calculate delay relative to control
      #self.VAFeven = np.linspace(45000.0, 209053.0, 2000)
      self.VAFeven = np.linspace(15000.0, 239053.0, 2000)
      self.timeOnVAFeven = np.interp(self.VAFeven, self.VAF[::-1], self.yrs[::-1])

      if not run['name'] in  ('ctrl', 'PIGLctrl'):
          steadyTimeOnVAFeven = runs[0]['data'].gs.timeOnVAFeven
          self.delay = self.timeOnVAFeven - steadyTimeOnVAFeven
          self.delay[np.nonzero(self.VAFeven > self.VAF.max())] = np.nan # remove nonsensical delays
          self.delay[np.nonzero(self.VAFeven < self.VAF.min())] = np.nan

class GIAoutputData:
   def __init__(self, run):
      # --------
      if run['name'] in  ('ctrl', 'PIGLctrl'):
         self.yrs = [0,]
         self.upliftVol = [0,]
      else:
         tic = time.perf_counter()
         f = netCDF4.Dataset(run['path'] + '/' + 'iceload_all.nc', 'r')
         nt = len(f.dimensions['Time'])
         self.yrs = np.arange(2015, 2015 + nt)
   
         x = f.variables['x'][:]
         y = f.variables['y'][:]
         dx = x[1] - x[0]
         [XX, YY] = np.meshgrid(x, y)
   
         bed = f.variables['bas']
         thk = f.variables['thk']

         toc = time.perf_counter()
         #print(f"GIA: Set up iceload_all.nc var objects in {toc - tic:0.4f} seconds")

   
         tic = time.perf_counter()
         self.upliftVol = np.zeros((nt,))
         self.upliftMean = np.zeros((nt,))
         self.upliftOceanVol = np.zeros((nt,))
         self.vaf = np.zeros((nt,))
   
         bed0 = bed[0,:,:]
         thk0 = thk[0,:,:]
   
         f2 = netCDF4.Dataset(run['path'] + '/' + 'uplift_GIA_all.nc', 'r')
         up = f2.variables['uplift']

         toc = time.perf_counter()
         #print(f"GIA: Set up uplift_GIA_all.nc var objects in {toc - tic:0.4f} seconds")

#         print(x.shape, y.shape, XX.shape, YY.shape, bed0.shape)
         for t in range(nt):
             tic = time.perf_counter()

             bedt = bed[t,:,:]
             thkt = thk[t,:,:]
             upt = up[t,:,:]
             #self.upliftVol[t] = ((bedt - bed0) * dx**2).sum()
             self.upliftVol[t] = (upt * dx**2).sum()
             self.upliftMean[t] = (upt).mean()
             

             hf = np.maximum(np.zeros(thkt.shape), -1.0 * bedt * rhow / rhoi)
             haf = np.maximum(np.zeros(thkt.shape), thkt - hf)
             self.vaf[t] = (dx**2 * haf).sum()

#             yrng = (YY>-520000.0)*(YY<-400000.0)
             [ind1, ind2] = np.where( (haf>1.0) * (haf<100.0) * (YY>-520000.0) * (YY<-400000.0))
#             if t % 100 == 0:
#                plt.matshow((YY>-520000.0) * (YY<-400000.0) * (haf>1.0) * (haf<100.0)); plt.show()
#             print(ind.shape, ind)
             mnGLx = XX[ind1, ind2].mean()
#             print("Time {}, mn GL x = {}".format(t, mnGLx))
             ocnMask = np.logical_or( (thkt<hf),
                                      (bed0 == 0.0) * (XX < mnGLx))  # outside TG domain and far enough grid-west
             #self.upliftOceanVol[t] = ((bedt - bed0) * dx**2 * ocnMask).sum()
             self.upliftOceanVol[t] = (upt * dx**2 * ocnMask).sum()

             toc = time.perf_counter()
             #print(f"GIA: processed year {t} in {toc - tic:0.4f} seconds")
         f.close()
         f2.close()

         if 'bad_uplift_ind' in run:
             # remove these
             self.yrs = np.delete(self.yrs, run['bad_uplift_ind'])
             self.upliftVol = np.delete(self.upliftVol, run['bad_uplift_ind'])
             self.upliftMean = np.delete(self.upliftMean, run['bad_uplift_ind'])
             self.vaf = np.delete(self.vaf, run['bad_uplift_ind'])
             self.upliftOceanVol = np.delete(self.upliftOceanVol, run['bad_uplift_ind'])

         if run['name'] == 'N1':
             # this run was kind of messed up
             self.yrs[134:] += 5


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

    axVAF.plot(yrs, VAF, label = run['legname'], color=run['color'])

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

# SL version
fig = plt.figure(33, facecolor='w')
axSL = fig.add_subplot(1, 1, 1)
for run in runs:
    print("Plotting run: " + run['name'])
    yrs = run['data'].gs.yrs
    VAF = run['data'].gs.VAF
    axSL.plot(yrs, GTtoSL(VAF-VAF[0]), label = run['legname'], color=run['color'])
axSL.legend(loc='best', ncol=1)
axSL.set_xlabel('Year')
axSL.set_ylabel('Sea-level rise (mm)')
axSL.grid(True)


# -------
# Bunch of stats
# -------
fig = plt.figure(2, facecolor='w', figsize=(8, 14))
nr = 6
axVAF = fig.add_subplot(nr, 1, 1)
#axGA = fig.add_subplot(nr, 1, 2)
axMelt = fig.add_subplot(nr, 1, 2)
#axFloatThk = fig.add_subplot(nr, 1, 4)
#axGrdThk = fig.add_subplot(nr, 1, 5)
axReduc = fig.add_subplot(nr, 1, 5)
axDelay = fig.add_subplot(nr, 1, 6)
axGLf = fig.add_subplot(nr, 1, 3)
axVAFrate = fig.add_subplot(nr, 1, 4)

for run in runs:
    print("Plotting run: " + run['name'])
    yrs = run['data'].gs.yrs
    VAF = run['data'].gs.VAF

    axVAF.plot(yrs, run['data'].gs.VAF, label = run['legname'], color=run['color'], linestyle=run['style'])
    #axGA.plot (yrs, run['data'].gs.GA, label = run['legname'], color=run['color'])
    axMelt.plot(yrs, run['data'].gs.melt, label = run['legname'], color=run['color'], linestyle=run['style'])
    #axMelt.plot(yrs, run['data'].gs.meltrate, label = run['legname'], color=run['color'])
    if not run['name'] in  ('ctrl', 'PIGLctrl'):
       axReduc.plot(run['data'].gs.resampYrs, run['data'].gs.reduction, label = run['legname'], color=run['color'], linestyle=run['style'])
       axDelay.plot(run['data'].gs.timeOnVAFeven, run['data'].gs.delay, label = run['legname'], color=run['color'], linestyle=run['style'])
    #axFloatThk.plot(yrs, run['data'].gs.floatThk, label = run['legname'], color=run['color'])
    #axGrdThk.plot(yrs, run['data'].gs.grdThk, label = run['legname'], color=run['color'])
    axGLf.plot(yrs, run['data'].gs.GLflux, label = run['legname'], color=run['color'], linestyle=run['style'])
    #axVAFrate.plot(yrs, run['data'].gs.VAFrate, label = run['legname'], color=run['color'])
    axVAFrate.plot(run['data'].gs.resampYrs[:-1], run['data'].gs.resampVAFrate, label = run['legname'], color=run['color'], linestyle=run['style'])

axVAF.legend(loc='best', ncol=1)
#axVAF.set_xlabel('Year')
axVAF.set_ylabel('VAF (Gt)')
axVAF.tick_params(bottom=True, top=True, left=True, right=False)
axVAF.grid(True)

axSLR=axVAF.twinx()
y1, y2=axVAF.get_ylim()
x1, x2=axVAF.get_xlim()
axSLR.set_ylim(GTtoSL(y1) - GTtoSL(VAF[0]), GTtoSL(y2) - GTtoSL(VAF[0]))
#axSLR.set_yticks( range(int(GTtoSL(y1)), int(GTtoSL(y2))) )
axSLR.set_ylabel('S.L. equiv. (mm)')
axSLR.set_xlim(x1, x2)
axSLR.tick_params(bottom=False, top=False, left=False, right=True)


##axGA.legend(loc='best', ncol=1)
#axGA.set_xlabel('Year')
#axGA.set_ylabel('Grounded area (km$^2$)')
#axGA.tick_params(bottom=True, top=True, left=True, right=True)

#axMelt.set_xlabel('Year')
axMelt.set_ylabel('Ice shelf melt\nrate (Gt yr${^-1}$)')
#axMelt.set_ylim((0.0, 2000.0))
axMelt.tick_params(bottom=True, top=True, left=True, right=True)
axMelt.grid(True)

#axReduc.set_xlabel('Year')
axReduc.set_ylabel('Reduction in VAF loss\nfrom control (%)')
axReduc.tick_params(bottom=True, top=True, left=True, right=True)
axReduc.grid(True)

axDelay.set_xlabel('Year')
axDelay.set_ylabel('Delay in VAF loss\nfrom control (yr)')
axDelay.tick_params(bottom=True, top=True, left=True, right=True)
axDelay.grid(True)

#axFloatThk.set_xlabel('Year')
#axFloatThk.set_ylabel('Mean floating ice thickness (m)')

#axGrdThk.set_xlabel('Year')
#axGrdThk.set_ylabel('Mean grounded ice thickness (m)')

#axGLf.set_xlabel('Year')
axGLf.set_ylabel('Grounding line\nflux (Gt yr${^-1}$)')
axGLf.tick_params(bottom=True, top=True, left=True, right=True)
axGLf.grid(True)

#axVAFrate.set_xlabel('Year')
axVAFrate.set_ylabel('VAF rate (Gt yr${^-1}$)')
axVAFrate.tick_params(bottom=True, top=True, left=True, right=True)
axVAFrate.grid(True)

# -------
# Bunch of stats vs GA
# -------
fig = plt.figure('vsGA', facecolor='w', figsize=(8, 14))
nr = 6
axVAF = fig.add_subplot(nr, 1, 1)
#axGA = fig.add_subplot(nr, 1, 2)
axMelt = fig.add_subplot(nr, 1, 2)
axFloatThk = fig.add_subplot(nr, 1, 5)
axGrdThk = fig.add_subplot(nr, 1, 4)
axGLf = fig.add_subplot(nr, 1, 3)

for run in runs:
    print("Plotting run: " + run['name'])
    GA = run['data'].gs.GAloss

    axVAF.plot(GA, run['data'].gs.VAF, label = run['legname'], color=run['color'])
    #axGA.plot (GA, run['data'].gs.GA, label = run['legname'], color=run['color'])
    axMelt.plot(GA, run['data'].gs.melt, label = run['legname'], color=run['color'])
    axFloatThk.plot(GA, run['data'].gs.floatThk, label = run['legname'], color=run['color'])
    axGrdThk.plot(GA, run['data'].gs.grdThk, label = run['legname'], color=run['color'])
    axGLf.plot(GA, run['data'].gs.GLflux, label = run['legname'], color=run['color'])

axVAF.legend(loc='best', ncol=1)
axVAF.set_xlabel('Grounded area loss (km$^2$)')
axVAF.set_ylabel('VAF (Gt)')
axVAF.tick_params(bottom=True, top=True, left=True, right=True)

axSLR=axVAF.twinx()
y1, y2=axVAF.get_ylim()
x1, x2=axVAF.get_xlim()
axSLR.set_ylim(GTtoSL(y1) - GTtoSL(VAF[0]), GTtoSL(y2) - GTtoSL(VAF[0]))
#axSLR.set_yticks( range(int(GTtoSL(y1)), int(GTtoSL(y2))) )
axSLR.set_ylabel('S.L. equiv. (mm)')
axSLR.set_xlim(x1, x2)
axSLR.tick_params(bottom=True, top=True, left=True, right=True)


#axGA.legend(loc='best', ncol=1)
#axGA.set_xlabel('Grounded area loss (km$^2$)')
#axGA.set_ylabel('Grounded area (km$^2$)')
#axGA.tick_params(bottom=True, top=True, left=True, right=True)

#axMelt.set_xlabel('Grounded area loss (km$^2$)')
axMelt.set_ylabel('Ice shelf melt\nrate (Gt yr${^-1}$)')
axMelt.set_ylim((0.0, 2000.0))
axMelt.tick_params(bottom=True, top=True, left=True, right=True)

axFloatThk.set_xlabel('Grounded area loss (km$^2$)')
axFloatThk.set_ylabel('Mean floating\nice thickness (m)')
axFloatThk.tick_params(bottom=True, top=True, left=True, right=True)

#axGrdThk.set_xlabel('Grounded area loss (km$^2$)')
axGrdThk.set_ylabel('Mean grounded ice\nthickness (m)')
axGrdThk.tick_params(bottom=True, top=True, left=True, right=True)

#axGLf.set_xlabel('Grounded area loss (km$^2$)')
axGLf.set_ylabel('Grounding line\nflux (Gt yr${^-1}$)')
axGLf.tick_params(bottom=True, top=True, left=True, right=True)

# --------
# GIA grid uplift time series figure
# --------
if doGIAfiles:
 fig = plt.figure('upGIA', facecolor='w', figsize=(5, 5))
 #axUp = fig.add_subplot(1, 2, 1)
 axUp2 = fig.add_subplot(1, 1, 1)
 for run in runs:
     if run['name'] in  ('ctrl', 'PIGLctrl'):
         continue
     print("Plotting run: " + run['name'])
     yrs = run['data'].GIA.yrs
     upliftMean = run['data'].GIA.upliftMean
     upliftVol = run['data'].GIA.upliftVol / (362.0e6 * 1000.0**2) * 1000.0 # mm
     upliftOceanVol = run['data'].GIA.upliftOceanVol / (362.0e6 * 1000.0**2) * 1000.0 # mm
     vaf = run['data'].GIA.vaf / (362.0e6 * 1000.0**2) * 1000.0 * rhoi/rhow # mm ocn
     bary = vaf[0] - vaf
 
     #axUp.plot(yrs, upliftVol, label = run['legname'], color=run['color'])
     #axUp.plot(yrs, upliftOceanVol, '--', color=run['color'])
 
     axUp2.plot(yrs, bary, '-', label = run['legname'], color=run['color'])
     axUp2.plot(yrs, upliftOceanVol, '--', color=run['color'])
     axUp2.fill_between(yrs, bary, bary + upliftOceanVol, color=run['color'], alpha=0.3)
 
 #axUp.legend(loc='best', ncol=1)
 #axUp.set_xlabel('Year')
 #axUp.set_ylabel('Uplift volume (mm SLE)')
 #axUp.tick_params(bottom=True, top=True, left=True, right=True)
 
 axUp2.legend(loc='best', ncol=1)
 axUp2.set_xlabel('Year')
 axUp2.set_ylabel('Sea level rise equivalent (mm)')
 axUp2.tick_params(bottom=True, top=True, left=True, right=True)
 axUp2.grid(True)

 fig9 = plt.figure('ocean error', facecolor='w', figsize=(5, 5))
 axUp9 = fig9.add_subplot(1, 1, 1)
 for run in runs:
     if run['name'] in  ('ctrl', 'PIGLctrl'):
         continue
     print("Plotting run: " + run['name'])
     yrs = run['data'].GIA.yrs
     upliftMean = run['data'].GIA.upliftMean
     upliftVol = run['data'].GIA.upliftVol / (362.0e6 * 1000.0**2) * 1000.0 # mm
     upliftOceanVol = run['data'].GIA.upliftOceanVol / (362.0e6 * 1000.0**2) * 1000.0 # mm
     vaf = run['data'].GIA.vaf / (362.0e6 * 1000.0**2) * 1000.0 * rhoi/rhow # mm ocn
     bary = vaf[0] - vaf
 
     #axUp.plot(yrs, upliftVol, label = run['legname'], color=run['color'])
     #axUp.plot(yrs, upliftOceanVol, '--', color=run['color'])
 
     axUp9.plot(yrs,upliftOceanVol/(bary), '-', label = run['legname'], color=run['color'])
 axUp9.legend(loc='best', ncol=1)
 axUp9.set_xlabel('Year')
 axUp9.set_ylabel('ocean disp/bary')
 axUp9.tick_params(bottom=True, top=True, left=True, right=True)
 axUp9.grid(True)





# --------
# uplift time series figure
# --------
if doUpliftTS:
 fig = plt.figure('up', facecolor='w')
 axUp = fig.add_subplot(1, 2, 1)
 axUp2 = fig.add_subplot(1, 2, 2)
 for run in runs:
     print("Plotting run: " + run['legname'])
     yrs = run['data'].out.yrs
     maxUplift = run['data'].out.maxUplift
     maxGrdUplift = run['data'].out.maxGrdUplift
     vaf = run['data'].out.vaf
 
     axUp.plot(yrs, maxUplift, label = run['legname'], color=run['color'])
     axUp.plot(yrs, maxGrdUplift, '--', color=run['color'])
 
     axUp2.plot(vaf[0]-vaf, maxUplift, label = run['legname'], color=run['color'])
     axUp2.plot(vaf[0]-vaf, maxGrdUplift, '--', color=run['color'])
 
 axUp.legend(loc='best', ncol=1)
 axUp.set_xlabel('Year')
 axUp.set_ylabel('Uplift (m)')
 
 axUp2.legend(loc='best', ncol=1)
 axUp2.set_xlabel('VAF loss (Gt)')
 axUp2.set_ylabel('Uplift (m)')


# --------
plt.show()
#plt.savefig('TG_time_series.png', dpi=300)

