from parcels import FieldSet, ParticleSet, AdvectionRK4, ErrorCode, plotTrajectoriesFile
from parcels import JITParticle, ScipyParticle, Variable

import stuckparticles_analysis as sta
import stuckparticles_class as stc
import stuckparticles_coasts as stco
import stuckparticles_field as stf
import stuckparticles_general as stg
import stuckparticles_grid as stgr
import stuckparticles_particles as stp
import stuckparticles_plot as stpl

# from multiprocessing.dummy import Pool as pool
from multiprocessing import Pool as pool

import os
import errno
import copy
import numpy as np
try:
    import matplotlib.pyplot as plt
except:
    pass

from datetime import timedelta


fieldlocation = "/data2/imau/oceanparcels/hydrodynamic_data/GlobCurrent/v2p0/total_hs/all_00hrs/"
fieldfiles = "/data2/imau/oceanparcels/hydrodynamic_data/GlobCurrent/v2p0/total_hs/all_00hrs/20020101"
# fieldlocation = "GlobCurrent/"
# fieldfiles = "GlobCurrent/"

runtime = timedelta(days=150)

# lons = np.arange(-175, 175.1, 2.5)
# lats = np.arange(-75, 75.1, 2.5)
lons = np.arange(-175, 175.1, 1)
lats = np.arange(-75, 75.1, 1)

time_stuck = 5. * 24*60*60
time_moving = 5. * 24*60*60


def timestepTest(timestep):
    print "Start timestepTest() for", timestep

    dt = timedelta(minutes=timestep)

    global lons
    global lats
    global runtime
    global fieldlocation
    global fieldfiles
    global time_stuck

    savefile = "timestep-{}.npz".format(timestep)

    ## FIELDS
    fset = stf.getFieldsetGlobCurrent(filelocation=fieldlocation)
    coasts = stco.findCoasts(filelocation=fieldfiles)

    stco.addCoasts(fset, coasts, constant=0)

    fset.add_periodic_halo(zonal=True)

    ## PARTICLES
    lons_mesh, lats_mesh = np.meshgrid(lons, lats)
    pset = ParticleSet(fieldset=fset, pclass=stc.stuckParticle, lon=lons_mesh, lat=lats_mesh)

    stf.removeLandParticles(particleset=pset, filename=fieldfiles)

    ## KERNELS
    kernels = AdvectionRK4 + pset.Kernel(stc.checkVelocity) + pset.Kernel(stg.periodicBC) + pset.Kernel(stco.returnFromCoast)

    pset.execute(
        kernels,
        runtime=runtime,
        dt=dt,
        recovery={ErrorCode.ErrorOutOfBounds: stg.deleteParticle}
    )

    ## RESULTS
    stuckparticles = 0
    for p in pset:
        if p.time_stuck > time_stuck:
            stuckparticles += 1

    sta.exportParticleData(fset, pset, savefile=savefile)
    print "Stop timestep", timestep, "with", stuckparticles
    return [timestep, stuckparticles]


def constantPush(constant):
    print "Start constantPush() for", constant

    global runtime
    global timestep_global
    global fieldlocation
    global fieldfiles
    global time_stuck

    dt = timedelta(minutes=timestep_global) #timedelta(minutes=10)
    savefile = "constantpush-{}-{}.npz".format(timestep_global, constant)

    ## FIELDS
    fset = stf.getFieldsetGlobCurrent(filelocation=fieldlocation)
    coasts = stco.findCoasts(filelocation=fieldfiles)

    stco.addCoasts(fset, coasts, constant=constant)

    fset.add_periodic_halo(zonal=True)

    ## PARTICLES
    lons_mesh, lats_mesh = np.meshgrid(lons, lats)
    pset = ParticleSet(fieldset=fset, pclass=stc.stuckParticle, lon=lons_mesh, lat=lats_mesh)

    stf.removeLandParticles(particleset=pset, filename=fieldfiles)

    ## KERNELS
    kernels = AdvectionRK4 + pset.Kernel(stc.checkVelocity) + pset.Kernel(stg.periodicBC) + pset.Kernel(stco.returnFromCoast)

    pset.execute(
        kernels,
        runtime=runtime,
        dt=dt,
        recovery={ErrorCode.ErrorOutOfBounds: stg.deleteParticle}
    )

    ## RESULTS
    stuckparticles = 0
    for p in pset:
        if p.time_stuck > time_stuck:
            stuckparticles += 1

    sta.exportParticleData(fset, pset, savefile=savefile)
    print "Stop constant", constant, "with", stuckparticles
    return [constant, stuckparticles]

pool_0 = pool()
timesteps_0 = np.arange(1, 20, 1, dtype=np.int)
results_0 = pool_0.map(timestepTest, timesteps_0)
np.savez_compressed("global-timestep-results", data=results_0)


constants_ = [0] + np.logspace(1, -7, 25).tolist() + [0]

timestep_global = 5
pool_5 = pool()
results_05 = pool_5.map(constantPush, constants_)
np.savez_compressed("global-constant-results-05", data=results_05)

timestep_global = 10
pool_10 = pool()
results_10 = pool_10.map(constantPush, constants_)
np.savez_compressed("global-constant-results-10", data=results_10)

timestep_global = 15
pool_15 = pool()
results_15 = pool_15.map(constantPush, constants_)
np.savez_compressed("global-constant-results-15", data=results_15)
