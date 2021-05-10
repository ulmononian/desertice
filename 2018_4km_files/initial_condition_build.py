#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 11:34:28 2019

@author: cbook
"""

import netCDF4


input = netCDF4.Dataset('thwaites.4km.cleaned.nc','r'); 
thk = input.variables['thickness'][:]
bas = input.variables['bedTopography'][:]

fout = netCDF4.Dataset('initial_condition.nc','w', format='NETCDF3_CLASSIC')

#thk = fout.createvariable('bedTopography')