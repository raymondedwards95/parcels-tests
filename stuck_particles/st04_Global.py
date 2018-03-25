""" Observing particles using kernels in stuckparticles_class.py """

from parcels import FieldSet, ParticleSet, Variable, JITParticle
from parcels import AdvectionRK4, plotTrajectoriesFile, ErrorCode, ScipyParticle

import stuckparticles as st
import stuckparticles_class as stc
from stuckparticles_class import DeleteParticle

import numpy as np
import math

from datetime import timedelta, datetime
from operator import attrgetter

import argparse
import os, errno


## Parameters
print "Start"
num_days = 75
subdomain = False
lon = np.linspace(-175, 175, num=71)
lat = np.linspace(-75, 75, num=31)

print "lon:", lon
print "lat:", lat

## Process parameters
# filename = "plots_04/p_{:05.3f}_{:05.3f}/".format(lon, lat)
# try:
#     os.makedirs(filename)
# except OSError as e:
#     if e.errno != errno.EEXIST:
#         raise

if subdomain:
    if np.size(lon) == 1:
        lon_0, lon_1 = int(np.floor(lon-5)), int(np.ceil(lon+5))
    else:
        lon_0, lon_1 = int(np.round(np.min(lon))), int(np.round(np.max(lon)))

    if np.size(lat) == 1:
        lat_0, lat_1 = int(np.floor(lat-5)), int(np.ceil(lat+5))
    else:
        lat_0, lat_1 = int(np.round(np.min(lat))), int(np.round(np.max(lat)))

    lon_range = range((lon_0-5+180)*4-1, (lon_1+5+180)*4+1)
    lat_range = range((lat_0-5+80)*4-1, (lat_1+5+80)*4+1)

    ind = {"lon": lon_range,
           "lat": lat_range}
else:
    ind = {}


## Get field
print "Get field"
filenames = {"U": "/data2/imau/oceanparcels/hydrodynamic_data/GlobCurrent/v2p0/total_hs/all_00hrs/20*.nc",
             "V": "/data2/imau/oceanparcels/hydrodynamic_data/GlobCurrent/v2p0/total_hs/all_00hrs/20*.nc"}
var = {"U": "eastward_eulerian_current_velocity",
       "V": "northward_eulerian_current_velocity"}
dim = {"lat": "lat",
       "lon": "lon",
       "time": "time"}
fieldset = FieldSet.from_netcdf(filenames, variables=var, dimensions=dim, indices=ind)


## Create particle
print "Create particles"
lon_mesh, lat_mesh = np.meshgrid(lon, lat)

pset = ParticleSet(fieldset=fieldset, pclass=stc.stuckParticle, lon=lon_mesh, lat=lat_mesh)
p0 = len(pset)

stc.removeLandParticles(fieldset, pset)
p1 = len(pset)

## Create kernel
print "Create additional kernel"
kernels = AdvectionRK4

k_velocity = pset.Kernel(stc.checkVelocity)
kernels += k_velocity

if subdomain is False:
    fieldset.add_periodic_halo(zonal=True)
    k_boundary = pset.Kernel(stc.periodicBC)
    kernels += k_boundary
else:
    def emptyKernel(particle, fieldset, time, dt):
        0.
    k_boundary = pset.Kernel(emptyKernel)


## Advect
print "Advecting with {} particles".format(len(pset))
pset.execute(
    AdvectionRK4 + k_velocity + k_boundary,
    runtime=timedelta(days=num_days),
    dt=timedelta(minutes=10),
    output_file=pset.ParticleFile(name="st04/trajectory_04a", outputdt=timedelta(hours=3)),
    recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle}
)
print "Advection finished"


## Save data
print "Export data"
data = stc.exportParticleData(fieldset, pset, velocities=True, savefile="st04_data_p{}_d{}".format(p1, num_days))


print "End"
