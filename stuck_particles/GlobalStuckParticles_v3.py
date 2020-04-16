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
# from multiprocessing import Pool as pool

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

runtime = timedelta(days=200)
# runtime = timedelta(days=20)

lons = np.arange(-175, 175.1, 5)
lats = np.arange(-65, 65.1, 5)
# lons = np.arange(-175, 175.1, 1)
# lats = np.arange(-75, 75.1, 1)
# lons = np.arange(-175, 175.1, 15)
# lats = np.arange(-65, 65.1, 15)

time_stuck = 5. * 24*60*60
time_moving = 5. * 24*60*60

output = False # save to pfile

def timestepTest(timestep):
    print "Start timestepTest() for", timestep

    dt = timedelta(minutes=timestep)

    # not needed:
    global lons
    global lats
    global runtime
    global fieldlocation
    global fieldfiles
    global time_stuck
    global time_moving
    global output

    savefile = "timestep-{}.npz".format(timestep)

    ## FIELDS
    fset = stf.getFieldsetGlobCurrent(filelocation=fieldlocation)
    coasts = stco.findCoasts(filelocation=fieldfiles)

    stco.addCoasts(fset, coasts, constant=0)

    fset.add_periodic_halo(zonal=True)

    ## PARTICLES
    lons_mesh, lats_mesh = np.meshgrid(lons, lats)
    pset = ParticleSet(fieldset=fset, pclass=stc.StuckCoastParticle, lon=lons_mesh, lat=lats_mesh)

    stf.removeLandParticles(particleset=pset, filename=fieldfiles)

    ## KERNELS
    kernels = AdvectionRK4 + pset.Kernel(stc.checkVelocity) + pset.Kernel(stg.periodicBC) + pset.Kernel(stco.returnFromCoast_A) + pset.Kernel(stc.particleOnCoast)

    if output:
        pset.execute(
            kernels,
            runtime=runtime,
            dt=dt,
            recovery={ErrorCode.ErrorOutOfBounds: stg.deleteParticle},
            output_file=pset.ParticleFile(name="pfile-"+str(timestep), outputdt=timedelta(hours=12))
        )
    else:
        pset.execute(
            kernels,
            runtime=runtime,
            dt=dt,
            recovery={ErrorCode.ErrorOutOfBounds: stg.deleteParticle}
        )


    ## RESULTS
    stuckparticles = 0
    coastparticles = 0
    oceanparticles = 0
    halfparticles_coast = 0
    halfparticles_ocean = 0

    total_time_coast = np.array([])
    total_time_ocean = np.array([])

    for p in pset:
        if p.time_stuck > time_stuck:
            stuckparticles += 1
        if p.number_on_coast == 0 or p.total_time_coast == 0.:
            oceanparticles += 1
        elif p.number_in_ocean == 0 or p.total_time_ocean == 0.:
            coastparticles += 1
        elif p.current_time_ocean > 0.:
            halfparticles_ocean += 1
        elif p.current_time_coast > 0.:
            halfparticles_coast += 1

        total_time_coast = np.append(total_time_coast, p.total_time_coast)
        total_time_ocean = np.append(total_time_ocean, p.total_time_ocean)


    mean_coast = np.mean(np.trim_zeros(np.sort(total_time_coast)))
    mean_coast_all = np.mean(np.sort(total_time_coast))
    mean_ocean = np.mean(np.trim_zeros(np.sort(total_time_ocean)))
    mean_ocean_all = np.mean(np.sort(total_time_ocean))


    sta.exportParticleData(fset, pset, savefile=savefile)

    print "Stop timestep", timestep, "with", stuckparticles

    return [timestep, stuckparticles, coastparticles, oceanparticles, halfparticles_coast, halfparticles_ocean, len(pset), mean_coast, mean_coast_all, mean_ocean, mean_ocean_all]


def constantPush(constant):
    print "Start constantPush() for", constant

    # not needed:
    global lons
    global lats
    global runtime
    global fieldlocation
    global fieldfiles
    global time_stuck
    global time_moving
    global output

    dt = timedelta(minutes=timestep_global) #timedelta(minutes=10)
    savefile = "constantpush-{}-{}.npz".format(timestep_global, constant)

    ## FIELDS
    fset = stf.getFieldsetGlobCurrent(filelocation=fieldlocation)
    coasts = stco.findCoasts(filelocation=fieldfiles)

    stco.addCoasts(fset, coasts, constant=constant)

    fset.add_periodic_halo(zonal=True)

    ## PARTICLES
    lons_mesh, lats_mesh = np.meshgrid(lons, lats)
    pset = ParticleSet(fieldset=fset, pclass=stc.StuckCoastParticle, lon=lons_mesh, lat=lats_mesh)

    stf.removeLandParticles(particleset=pset, filename=fieldfiles)

    ## KERNELS
    kernels = AdvectionRK4 + pset.Kernel(stc.checkVelocity) + pset.Kernel(stg.periodicBC) + pset.Kernel(stco.returnFromCoast_A) + pset.Kernel(stc.particleOnCoast)

    if output:
        pset.execute(
            kernels,
            runtime=runtime,
            dt=dt,
            recovery={ErrorCode.ErrorOutOfBounds: stg.deleteParticle},
            output_file=pset.ParticleFile(name="pfile-"+str(timestep)+"-"+str(constant), outputdt=timedelta(hours=12))
        )
    else:
        pset.execute(
            kernels,
            runtime=runtime,
            dt=dt,
            recovery={ErrorCode.ErrorOutOfBounds: stg.deleteParticle}
        )

    ## RESULTS
    stuckparticles = 0
    coastparticles = 0
    oceanparticles = 0
    halfparticles_coast = 0
    halfparticles_ocean = 0

    total_time_coast = np.array([])
    total_time_ocean = np.array([])

    for p in pset:
        if p.time_stuck > time_stuck:
            stuckparticles += 1
        if p.number_on_coast == 0 or p.total_time_coast == 0.:
            oceanparticles += 1
        elif p.number_in_ocean == 0 or p.total_time_ocean == 0.:
            coastparticles += 1
        elif p.current_time_ocean > 0.:
            halfparticles_ocean += 1
        elif p.current_time_coast > 0.:
            halfparticles_coast += 1

        total_time_coast = np.append(total_time_coast, p.total_time_coast)
        total_time_ocean = np.append(total_time_ocean, p.total_time_ocean)


    mean_coast = np.mean(np.trim_zeros(np.sort(total_time_coast)))
    mean_coast_all = np.mean(np.sort(total_time_coast))
    mean_ocean = np.mean(np.trim_zeros(np.sort(total_time_ocean)))
    mean_ocean_all = np.mean(np.sort(total_time_ocean))


    sta.exportParticleData(fset, pset, savefile=savefile)

    print "Stop constant", constant, "with", stuckparticles

    return [constant, stuckparticles, coastparticles, oceanparticles, halfparticles_coast, halfparticles_ocean, len(pset), mean_coast, mean_coast_all, mean_ocean, mean_ocean_all]

# timestepTest(10)

# pool_0 = pool()

timesteps_a = np.arange(1, 30.1, 1)
timesteps_b = np.arange(5, 60.1, 5)/60.
timesteps_c = np.arange(30, 90.1, 5)
timesteps_0 = np.round(np.sort(np.unique(np.concatenate((timesteps_a, timesteps_b, timesteps_c)))), 3)

# # results_a = pool_0.map(timestepTest, timesteps_a)
# results_a = map(timestepTest, timesteps_a)
# results_a = np.array(list(results_a))
# np.savez_compressed("global-timestep-results-a", data=results_a)
# # results_b = pool_0.map(timestepTest, timesteps_b)
# results_b = map(timestepTest, timesteps_b)
# results_b = np.array(list(results_b))
# np.savez_compressed("global-timestep-results-b", data=results_b)
# # results_c = pool_0.map(timestepTest, timesteps_c)
# results_c = map(timestepTest, timesteps_c)
# results_c = np.array(list(results_c))
# np.savez_compressed("global-timestep-results-c", data=results_c)
# results_0 = pool_0.map(timestepTest, timesteps_0)
results_0 = map(timestepTest, timesteps_0)
results_0 = np.array(list(results_0))
np.savez_compressed("global-timestep-results", data=results_0)


constants_ = [0] + np.logspace(3, -3, 25).tolist() + [0]

timestep_global = 5
# pool_5 = pool()
# results_05 = pool_5.map(constantPush, constants_)
results_05 = map(constantPush, constants_)
results_05 = np.array(list(results_05))
np.savez_compressed("global-constant-results-05", data=results_05)

timestep_global = 10
# pool_10 = pool()
# results_10 = pool_10.map(constantPush, constants_)
results_10 = map(constantPush, constants_)
results_10 = np.array(list(results_10))
np.savez_compressed("global-constant-results-10", data=results_10)

timestep_global = 15
# pool_15 = pool()
# results_15 = pool_15.map(constantPush, constants_)
results_15 = map(constantPush, constants_)
results_15 = np.array(list(results_15))
np.savez_compressed("global-constant-results-15", data=results_15)


output = True
print " ===== END ===== "
