from parcels import FieldSet, ParticleSet, Variable, JITParticle
from parcels import AdvectionRK4, plotTrajectoriesFile, ErrorCode, ScipyParticle

import numpy as np
import math

from datetime import timedelta, datetime
from operator import attrgetter

def DeleteParticle(particle, fieldset, time, dt):
    p = particle
    p.delete()
    print "ErrorOutOfBounds --> Delete Particle {} at ({}, {})".format(p.id, p.lon, p.lat)

filenames = {"U": "GlobCurrent_Agulhas/20*.nc",
             "V": "GlobCurrent_Agulhas/20*.nc"}
variables = {"U": "eastward_eulerian_current_velocity",
             "V": "northward_eulerian_current_velocity"}
dimensions = {"lat": "lat",
              "lon": "lon",
              "time": "time"}
fieldset = FieldSet.from_netcdf(filenames, variables, dimensions)

lon, lat = [18, 18], [-32.5, -33]

min_lon, max_lon = min(lon)-1, max(lon)+1
min_lat, max_lat = min(lat)-1, max(lat)+1
dom = [max_lat+1, min_lat-1, max_lon+1, min_lon-1]

pset = ParticleSet(fieldset=fieldset, pclass=JITParticle, lon=lon, lat=lat)

number = 50
for num in range(number):
    # First plot the particles
    print pset
    pset.show(savefile="01_plots/particles"+str(num).zfill(3), field="vector", domain=dom, vmax=2.0)

    # Then advect the particles for 6 hours
    pset.execute(AdvectionRK4,
                 runtime=timedelta(hours=24),
                 dt=timedelta(minutes=10),
                 recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

print pset
pset.show(savefile="01_plots/particles"+str(number).zfill(3), field="vector", domain=dom, vmax=2.0)
