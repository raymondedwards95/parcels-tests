""" Simple script to test with GlobCurrent
The files *.nc are in a folder GlobCurrent/
Resolution is 1/4
Range is [80, -80] (lat) and [-180, 180] (lon)
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

exportname = "globcur_01"


### Import field
filename = "201402*.nc"
location = "GlobCurrent/"

filenames = {"U": location+filename,
             "V": location+filename}
variables = {"U": "eastward_eulerian_current_velocity",
             "V": "northward_eulerian_current_velocity"}
dimensions = {"lat": "lat",
              "lon": "lon",
              "time": "time"}
indices = {}
fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, indices, allow_time_extrapolation=True)
if test: print "fieldset created"


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
                 dt=timedelta(minutes=30)
                 )

    if plot:
        pset.show(field=fieldset.U,
                  particles=False,
                  savefile=exportname+"_0U"
                  )
        pset.show(field=fieldset.V,
                  particles=False,
                  savefile=exportname+"_0V"
                  )
        pset.show(field='vector',
                  land=True,
                  particles=False,
                  vmax=3,
                  savefile=exportname+"_1"
                  )

    print pset
    print


if test:
    print "End of script"
    print
