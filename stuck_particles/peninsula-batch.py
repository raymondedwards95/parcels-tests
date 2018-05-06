## FROM NOTEBOOK Peninsula stuck particles timestep test
from parcels import FieldSet, ParticleSet, AdvectionRK4, ErrorCode, plotTrajectoriesFile
from parcels import JITParticle

import stuckparticles_analysis as sta
import stuckparticles_class as stc
import stuckparticles_coasts as stco
import stuckparticles_field as stf
import stuckparticles_general as stg
import stuckparticles_grid as stgr
import stuckparticles_particles as stp
import stuckparticles_plot as stpl

from peninsula_fieldset import peninsula_fieldset

import os
import errno
import copy
import numpy as np

from datetime import timedelta


def peninsula_timestep_test(timestep):
    print "Start timestep", timestep

    xdim, ydim = 50, 50
    stucktime = 0
    runtime = timedelta(days=2)

    p_lons = [0.2]
    p_lats = np.linspace(0.01, 0.2, num=201)
    p_lon_mesh, p_lat_mesh = np.meshgrid(p_lons, p_lats)

    fset = peninsula_fieldset(xdim, ydim, mesh='spherical')
    fset.add_periodic_halo(zonal=True)
    pset = ParticleSet(fieldset=fset, pclass=stc.stuckParticle, lon=p_lon_mesh, lat=p_lat_mesh)

    stf.removeLandParticles(fieldset=fset, particleset=pset)

    kernels = AdvectionRK4 + pset.Kernel(stc.checkVelocity) + pset.Kernel(stg.periodicBC)

    pset.execute(
        kernels,
        runtime=runtime,
        dt=timedelta(minutes=timestep),
        recovery={ErrorCode.ErrorOutOfBounds: stg.deleteParticle}
    )

    stuckparticles = 0
    for p in pset:
        if p.time_stuck > stucktime:
            stuckparticles += 1

    print "Stop timestep", timestep, "with", stuckparticles
    return [timestep, stuckparticles]


def peninsula_timestep_test_coasts(timestep):
    print "Start timestep", timestep

    xdim, ydim = 50, 50
    stucktime = 0
    runtime = timedelta(days=2)

    p_lons = [0.2]
    p_lats = np.linspace(0.01, 0.2, num=201)
    p_lon_mesh, p_lat_mesh = np.meshgrid(p_lons, p_lats)

    fset = peninsula_fieldset(xdim, ydim, mesh='spherical')
    coasts = stco.findCoasts(fieldset=fset)
    stco.addCoasts(fieldset=fset, coastfields=coasts)

    fset.add_periodic_halo(zonal=True)
    pset = ParticleSet(fieldset=fset, pclass=stc.stuckParticle, lon=p_lon_mesh, lat=p_lat_mesh)

    stf.removeLandParticles(fieldset=fset, particleset=pset)

    kernels = AdvectionRK4 + pset.Kernel(stc.checkVelocity) + pset.Kernel(stg.periodicBC) + pset.Kernel(stco.returnFromCoast)

    pset.execute(
        kernels,
        runtime=runtime,
        dt=timedelta(minutes=timestep),
        recovery={ErrorCode.ErrorOutOfBounds: stg.deleteParticle}
    )

    stuckparticles = 0
    for p in pset:
        if p.time_stuck > stucktime:
            stuckparticles += 1

    print "Stop timestep", timestep, "with", stuckparticles
    return [timestep, stuckparticles]


## Normal
timesteps_0 = np.arange(10, 1201, 10, dtype=np.int)
results_0 = map(peninsula_timestep_test, timesteps_0)

np.savez_compressed("peninsula-timestep-results", data=results_0)


## Coasts
timesteps_0 = np.arange(10, 1201, 10, dtype=np.int)
results_0 = map(peninsula_timestep_test_coasts, timesteps_0)

np.savez_compressed("peninsula-timestep-results-coasts", data=results_0)
