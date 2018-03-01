""" Test indices with OSCAR
The files oscar_vel*.nc are in a folder OSCAR/
Resolution is 1/3
Default range is [80, -80] (lat) and [20, 420] (lon)
"""
### Import modules, classes and functions
from parcels import FieldSet, ParticleSet
from parcels import AdvectionRK4
from parcels import JITParticle, ScipyParticle

import numpy as np
import math

from datetime import timedelta, datetime


### Tests
test = True
plot = True
particle = True
lists = False

exportname = "oscar_02"


### Find indices
min_lon, max_lon = 110, 170 # default 20, 420
min_lat, max_lat = -20, 50 # default -80, 80

lon_range = range((min_lon-20)*3, (max_lon-20)*3+1)
lat_range = range((-max_lat+80)*3, (-min_lat+80)*3+1)


### Import field
filename = "oscar_vel*.nc"
location = "OSCAR/"

filenames = {"U": location+filename,
             "V": location+filename}
variables = {"U": "u",
             "V": "v"}
dimensions = {"lat": "latitude",
              "lon": "longitude",
              "time": "time",
              "depth": "depth"}
indices = {'lon': lon_range,
           'lat': lat_range}
fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, indices, allow_time_extrapolation=True)
if test: print "fieldset created"


### Check lons and lats
if list:
    print "Parts of fieldset.*.data:"
    print fieldset.U.data[0,0,120:125]
    print fieldset.V.data[0,0,120:125]
    print
    print fieldset.U.data[0,0,-125:-120]
    print fieldset.V.data[0,0,-125:-120]
    print
    print "fieldset.U.grid.lat.shape", fieldset.U.grid.lat.shape
    print "fieldset.U.grid.lon.shape", fieldset.U.grid.lon.shape
    print "fieldset.U.data.shape    ", fieldset.U.data.shape
    print "fieldset.V.grid.lat.shape", fieldset.V.grid.lat.shape
    print "fieldset.V.grid.lon.shape", fieldset.V.grid.lon.shape
    print "fieldset.V.data.shape    ", fieldset.V.data.shape

if test:
    print
    print "fieldset.U.lon", fieldset.U.lon[:3], "...", fieldset.U.lon[-3:]
    print "fieldset.U.grid.lon", fieldset.U.grid.lon[:3], "...", fieldset.U.grid.lon[-3:]
    print
    print "fieldset.U.lat", fieldset.U.lat[:3], "...", fieldset.U.lat[-3:]
    print "fieldset.U.grid.lat", fieldset.U.grid.lat[:3], "...", fieldset.U.grid.lat[-3:]
    print
    print "fieldset.V.lon", fieldset.V.lon[:3], "...", fieldset.V.lon[-3:]
    print "fieldset.V.grid.lon", fieldset.V.grid.lon[:3], "...", fieldset.V.grid.lon[-3:]
    print
    print "fieldset.V.lat", fieldset.V.lat[:3], "...", fieldset.V.lat[-3:]
    print "fieldset.V.grid.lat", fieldset.V.grid.lat[:3], "...", fieldset.V.grid.lat[-3:]
    print
    print "lons should be ascending, from {} to {}".format(min_lon, max_lon)
    print "lats should be ascending, from {} to {}".format(min_lat, max_lat)
    print


### Create plots
if particle:
    pset = ParticleSet(fieldset=fieldset,
                       pclass=ScipyParticle,
                       lon=140, lat=0
                       )

    print pset
    print

    pset.execute(AdvectionRK4,
                 runtime=timedelta(days=1),
                 dt=timedelta(minutes=30),
                 # recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle}
                 )

    if plot:
        pset.show(field=fieldset.U,
                  land=True,
                  particles=False,
                  #domain=dom,
                  savefile=exportname+"_0"
                  )
        pset.show(field='vector',
                  land=True,
                  particles=False,
                  #domain=dom,
                  vmax=3,
                  savefile=exportname+"_1"
                  )

    print pset
    print


if test:
    print "End of script"
    print
