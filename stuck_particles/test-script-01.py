## Imports
from parcels import ParticleSet, AdvectionRK4, ErrorCode

import stuckparticles as st
import stuckparticles_class as stc
import stuckparticles_analysis as sta
import stuckparticles_field as stf

from stuckparticles_class import DeleteParticle

import matplotlib.pyplot as plt
import numpy as np
import math

import os
import errno
import copy

from datetime import timedelta


## Parameters
p_lon = np.linspace(-175, 175, num=71)
p_lat = np.linspace(-75, 75, num=31)

filename = "Particleset - 02/"

simulation_time = 75

# print 71*31

try:
    os.makedirs(filename)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise


## Fieldset
fieldset_0 = stf.getFieldsetGlobCurrent(filelocation="GlobCurrent/", time_extrapolation=True)
fieldset_1 = copy.deepcopy(fieldset_0)


## Find coasts
coast_fields_1 = stf.createCoastVelocities(fieldset=fieldset_1, factor=-1, abs=False, constant=0)
stf.exportCoastVelocities(coast_fields_1[0], coast_fields_1[1], filename="GlobCurrentCoasts_1")

stf.showCoast(coast_fields_1, type=np.bool)
stf.showCoast(coast_fields_1, type=np.float32)


## Add coasts to fieldset
fieldset_1 = stf.addGlobCurrentCoast(fieldset_1, coast_fields_1)


## Create particles
p_lon_mesh, p_lat_mesh = np.meshgrid(p_lon, p_lat)

pset_0 = ParticleSet(fieldset=fieldset_0, pclass=stc.stuckParticle, lon=p_lon_mesh, lat=p_lat_mesh)
pset_1 = ParticleSet(fieldset=fieldset_1, pclass=stc.stuckParticle, lon=p_lon_mesh, lat=p_lat_mesh)

## Remove landparticles
stc.removeLandParticles(fieldset_0, pset_0)
stc.removeLandParticles(fieldset_0, pset_1) # use same fieldset to remove the same particles

## Show results
print "\nNumber of particles in both psets:", len(pset_0), len(pset_1)


## Define kernels
kernels_0 = AdvectionRK4

kernels_1 = AdvectionRK4


## velocity
k_velocity_0 = pset_0.Kernel(stc.checkVelocity)
kernels_0 += k_velocity_0

k_velocity_1 = pset_1.Kernel(stc.checkVelocity)
kernels_1 += k_velocity_1


## periodic boundary conditions
fieldset_0.add_periodic_halo(zonal=True)
k_boundary_0 = pset_0.Kernel(stc.periodicBC)
kernels_0 += k_boundary_0

fieldset_1.add_periodic_halo(zonal=True)
k_boundary_1 = pset_1.Kernel(stc.periodicBC)
kernels_1 += k_boundary_1


## Advect
pset_0.execute(
    kernels_0,
    runtime=timedelta(days=simulation_time),
    dt=timedelta(minutes=10),
    output_file=pset_0.ParticleFile(name=filename+"trajectory_0", outputdt=timedelta(hours=3)),
    recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle}
)
print "Advection 0 finished"

pset_1.execute(
    kernels_1,
    runtime=timedelta(days=simulation_time),
    dt=timedelta(minutes=10),
    output_file=pset_1.ParticleFile(name=filename+"trajectory_1", outputdt=timedelta(hours=3)),
    recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle}
)
print "Advection 1 finished"


## Interpret results
stuck_particles_0 = 0
stuck_times_0 = []
free_particles_0 = 0

stuck_particles_1 = 0
stuck_times_1 = []
free_particles_1 = 0

for particle in pset_0:
    stuck_times_0.append(particle.time_stuck)
    if particle.time_stuck > 5*24*60*60:
        stuck_particles_0 += 1
    else:
        free_particles_0 += 1

for particle in pset_1:
    stuck_times_1.append(particle.time_stuck)
    if particle.time_stuck > 5*24*60*60:
        stuck_particles_1 += 1
    else:
        free_particles_1 += 1

print "Number of stuck particles:", stuck_particles_0, stuck_particles_1
print "Number of free particles: ", free_particles_0, free_particles_1

print "Total:", len(pset_0), len(pset_1)


## Histogram
print "total:", len(pset_0), len(pset_1)
print "stuck:", stuck_particles_0, stuck_particles_1

y_0 = np.array(sorted(stuck_times_0, reverse=True))/(24*60*60)
y_red_0 = np.trim_zeros(y_0)

y_1 = np.array(sorted(stuck_times_1, reverse=True))/(24*60*60)
y_red_1 = np.trim_zeros(y_1)

n_bins_0 = int(np.round(simulation_time/5.))

plt.figure()
plt.hist(y_red_0, n_bins_0, facecolor="green", alpha=0.5, label="pset_0")
plt.hist(y_red_1, n_bins_0, facecolor="red", alpha=0.5, label="pset_1")
plt.xlabel("days stuck")
plt.ylabel("number of particles")
plt.title("")
plt.legend()
plt.grid()
plt.show()


## Save data
data_0 = sta.exportParticleData(fieldset_0, pset_0, velocities=True, savefile=filename+"data_0")
data_1 = sta.exportParticleData(fieldset_1, pset_1, velocities=True, savefile=filename+"data_1")


## Extract stuck particles
stuck_0 = sta.extractStuckParticles(data_0, time_stuck=5, time_moving=0)
stuck_1 = sta.extractStuckParticles(data_1, time_stuck=5, time_moving=0)


## Print locations
sta.printLocations(stuck_0, initial=True)
sta.printLocations(stuck_1, initial=True)


## Plot locations
print stuck_particles_0
sta.plotLocations(stuck_0, coastfield=coast_fields_1)

print stuck_particles_1
sta.plotLocations(stuck_1, coastfield=coast_fields_1)


## Print grid velocities
sta.printGridVelocity(stuck_0, flux=True)
sta.printGridVelocity(stuck_1, flux=True)
