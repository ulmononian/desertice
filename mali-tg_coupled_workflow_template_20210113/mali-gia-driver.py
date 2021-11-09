#!/usr/bin/env python
"""
Author: Matt Hoffman
Date: February 28, 2019
Driver script to run gia code for data coupling with MALI
"""
import sys, os
import netCDF4
import numpy as np
import argparse
import configparser
# add path to giascript.py
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import giascript

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.description = __doc__
parser.add_argument("-r", dest="restart", action='store_true', help="Indicates the GIA model is to be restarted from a previous run.")
parser.add_argument("-i", "--input", dest="inputFile", help="name of ice load time-series input file", default="iceload.nc", metavar="FILENAME")
parser.add_argument("-t", "--timestep", dest="timeStep", help="timestep to be used by GIA model in years", default="0.1", type=float, metavar="DT")
parser.add_argument("-d", "--Duration", dest="Duration", help="run duration to be used by GIA model in years", default="1.0", type=float, metavar="T")
parser.add_argument("-ic", "--initialCondition", dest="initialConditionFile", help="name of thk, bas initial condition file", default="initial_condition.nc", metavar="FILENAME")

options = parser.parse_args()

maliRestart=options.restart

if maliRestart:
    print("this is a restart from a previous GIA run, so loading needed restart info.")
    fr = netCDF4.Dataset('gia_restart_data.nc','r')
    # TODO: check size matches size used below
    Uhatn_restart = fr.variables['Uhatn'][:]
    taf0hat_restart = fr.variables['taf0hat'][:]
    dLold_restart = fr.variables['dLold'][:]
    thk_Start = fr.variables['thk_Start'][:]
    bas_Start = fr.variables['bas_Start'][:]
    fr.close()
else:
    Uhatn_restart = None
    taf0hat_restart = None
    dLold_restart = None
    ic = netCDF4.Dataset(options.initialConditionFile, 'r')
    thk_Start = ic.variables['thk'][:]
    bas_Start = ic.variables['bas'][:]  # open special init cond file for original thk_start and bas_start
    ic.close()

# Open iceload file interpolated from MALI to GIA grid and get some grid info
f = netCDF4.Dataset(options.inputFile,'r')
x_data = f.variables['x'][:]
y_data = f.variables['y'][:]
Nx = len(x_data)
Ny = len(y_data)
xi, yj = np.meshgrid(np.arange(Nx), np.arange(Ny))
thk_End = f.variables['thk'][0,:,:]
bas_End = f.variables['bas'][0,:,:]
# nt = len(f.dimensions['Time'])
# note: Eventually may want to make the timekeeping more robust.
f.close()
T = options.Duration
dt = options.timeStep
nt = int(T / dt)
print("Using nt = {}".format(nt))
assert nt*dt==T, "dt does not divide evenly into duration."

print("Using dt = {} years".format(dt))

config = configparser.ConfigParser()
# Read config file
config.read('rheo.ini')

ekwargs = {}
if config.has_option('ekwargs', 'u2'):
   ekwargs['u2'] = float(config['ekwargs']['u2'])
if config.has_option('ekwargs', 'u1'):
   ekwargs['u1'] = float(config['ekwargs']['u1'])
if config.has_option('ekwargs', 'u'):
   ekwargs['u'] = float(config['ekwargs']['u'])
if config.has_option('ekwargs', 'h'):
   ekwargs['h']  = float(config['ekwargs']['h'])
ekwargs['D']  = float(config['ekwargs']['D'])

if ('u' in ekwargs and 'u2' in ekwargs) or ('u' in ekwargs and 'u1' in ekwargs):
   sys.exit("Error: If you specify u, you cannot specify u1 or u2")

if ('u1' in ekwargs and not 'u2' in ekwargs) or ('u2' in ekwargs and not 'u1' in ekwargs):
   sys.exit("Error: If you specify one of u1 and u2, you need to specify the other")

if 'u' in ekwargs and 'h' in ekwargs:
   sys.exit("Error: If you use a single layer viscosity (i.e. specify u instead of u1,u2), you should not include h.")

print("Using the following rheology parameters:", ekwargs)

maliForcing = giascript.maliForcing(thk_Start, bas_Start, thk_End, bas_End, nt)
buelerflux = giascript.BuelerTopgFlux(x_data, y_data, './', options.inputFile, 'blah', nt, dt, ekwargs, fac=3, read='netcdf_read', U0=Uhatn_restart, taf0=taf0hat_restart, dLold=dLold_restart, include_viscous=True, include_elastic=True, include_ocean=True, maliForcing=maliForcing)

# create a new GIA output file
fout = netCDF4.Dataset("uplift_GIA.nc", "w")
fout.createDimension('x', Nx)
fout.createDimension('y', Ny)
fout.createDimension('Time', size=None) # make unlimited dimension
xout = fout.createVariable('x', 'f', ('x',))
xout[:] = x_data
yout = fout.createVariable('y', 'f', ('y',))
yout[:] = y_data
tout = fout.createVariable('Time', 'f', ('Time',))
tout.units='year'
up_out = fout.createVariable('uplift', 'f', ('Time', 'y','x'))
uplift_rate = fout.createVariable('uplift_rate', 'f', ('Time', 'y','x'))


for i in range(nt):
    print("Starting time step {}".format(i))
    buelerflux._update_Udot(i)  # this actually runs the GIA model, using the MALI data from the iceload.nc file

up_out[0,:,:] = buelerflux.ifft2andcrop(buelerflux.Uhatn)
uplift_rate[0,:,:] = buelerflux.Udot
tout[0] = T

fout.close()

# Always write a restart file, clobbering previous one
fr = netCDF4.Dataset('gia_restart_data.nc','w')
fr.createDimension('x', Nx)
fr.createDimension('y', Ny)
xout = fr.createVariable('x', 'f', ('x',))
xout[:] = x_data
yout = fr.createVariable('y', 'f', ('y',))
yout[:] = y_data
Uhatn_restartVar = fout.createVariable('Uhatn', 'f', ('y','x'))
Uhatn_restartVar[:] = buelerflux.ifft2andcrop(buelerflux.Uhatn)
taf0hat_restartVar = fout.createVariable('taf0hat', 'f', ('y','x'))
taf0hat_restartVar[:] = buelerflux.ifft2andcrop(buelerflux.taf0hat)
dLold_restartVar = fout.createVariable('dLold', 'f', ('y','x'))
dLold_restartVar[:] = buelerflux.ifft2andcrop(buelerflux.dLhatold)
thk_Start = fout.createVariable('thk_Start', 'f', ('y','x'))
thk_Start[:] = thk_End
bas_Start = fout.createVariable('bas_Start', 'f', ('y','x'))
bas_Start[:] = bas_End

fr.close()

