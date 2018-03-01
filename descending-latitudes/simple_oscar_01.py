""" Simple script to test with OSCAR
The files oscar_vel*.nc are in a folder OSCAR/
Resolution is 1/3
Range is [80, -80] (lat) and [20, 420] (lon)
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

exportname = "oscar_01"


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
indices = {}
fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, indices, allow_time_extrapolation=True)
if test: print "fieldset created"


### Check lons and lats
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
    print "lons should be ascending"
    print "lats should be ascending"


### Create plots
if particle:
    pset = ParticleSet(fieldset=fieldset,
                       pclass=ScipyParticle,
                       lon=140, lat=0
                       )

    print pset

    pset.execute(AdvectionRK4,
                 runtime=timedelta(days=1),
                 dt=timedelta(minutes=30)
                 )

    if plot:
        pset.show(field=fieldset.U,
                  land=True,
                  particles=False,
                  savefile=exportname+"_0"
                  )
        pset.show(field='vector',
                  land=True,
                  particles=False,
                  domain=[80, -80, 400, 40],
                  vmax=3,
                  savefile=exportname+"_1"
                  )

    print pset


if test:
    print
    print "End of script"
