"""
##### IMPORT #####
"""
DEBUG = False # not working for whole file
if DEBUG: print "DEBUG: DEBUG IS ENABLED"

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

# from multiprocessing.dummy import Pool
from multiprocessing import Pool

import os
import errno
import copy
import numpy as np
try:
    import matplotlib.pyplot as plt
except:
    pass

from datetime import timedelta


"""
##### GLOBAL PARAMETERS #####
"""

### FIELDS + TIME-PARAMETERS
try:
    if DEBUG: print "DEBUG: try gemini"
    fieldlocation = "/data2/imau/oceanparcels/hydrodynamic_data/GlobCurrent/v2p0/total_hs/all_00hrs/"
    fieldfiles = "/data2/imau/oceanparcels/hydrodynamic_data/GlobCurrent/v2p0/total_hs/all_00hrs/20020101"

    fieldset = stf.getFieldsetGlobCurrent(filelocation=fieldlocation)
    coasts = stco.findCoasts(filelocation=fieldfiles)
    stco.exportCoasts(coasts, "coasts-globcurrent")

    runtime = timedelta(days=200)

    lons = np.arange(-175, 175.1, 2.5)
    lats = np.arange(-65, 65.1, 2.5)

    num = 3

    print "calculations are on Gemini"

except:
    if DEBUG: print "DEBUG: do local"
    fieldlocation = "GlobCurrent/"
    fieldfiles = "GlobCurrent/"

    fieldset = stf.getFieldsetGlobCurrent(filelocation=fieldlocation)
    fieldset.add_periodic_halo(zonal=True)
    coasts = stco.findCoasts(filelocation=fieldfiles)
    stco.exportCoasts(coasts, "coasts-globcurrent")

    runtime = timedelta(days=200) # 25?

    lons = np.arange(-175, 175.1, 2.5)
    lats = np.arange(-65, 65.1, 2.5)

    num = 7 # 7

    print "calculations are local (and short)"


### PARTICLES
# lons = np.arange(-175, 175.1, 5)
# lats = np.arange(-65, 65.1, 5)
# lons = np.arange(-175, 175.1, 15)
# lats = np.arange(-65, 65.1, 15)

lons_mesh, lats_mesh = np.meshgrid(lons, lats)

time_stuck = 5. * 24*60*60
time_moving = 5. * 24*60*60

output = False # save to pfile


"""
##### FUNCTION TIMESTEP #####
"""
def test_timestep(timestep):
    ### START
    print "\nStart test_timestep() for dt={}".format(timestep)

    dt = timedelta(minutes=timestep)
    savefile = "timestep-{:.3f}".format(timestep)

    ### ADJUST FIELD
    if DEBUG: print "DEBUG: get fset"
    fset = copy.copy(fieldset)
    coasts = stco.importCoasts("coasts-globcurrent.npz")
    if DEBUG: print "DEBUG: add coasts to fset"
    stco.addCoasts(fset, coasts, constant=0)

    ### PARTICLES
    pset = ParticleSet(fieldset=fset, pclass=stc.StuckCoastParticle, lon=lons_mesh, lat=lats_mesh)

    stf.removeLandParticles(particleset=pset, filename=fieldfiles)

    ### KERNELS
    kernels = AdvectionRK4 + pset.Kernel(stc.checkVelocity) + pset.Kernel(stg.periodicBC) + pset.Kernel(stco.returnFromCoast_A) + pset.Kernel(stc.particleOnCoast)

    ### RUN
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

    ### RESULTS
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

    ### EXPORT DATA
    sta.exportParticleData(fset, pset, savefile=savefile+"-data")

    print "Stop test_timestep for {} with {}, {}".format(timestep, stuckparticles, coastparticles)

    return [timestep, stuckparticles, coastparticles, oceanparticles, halfparticles_coast, halfparticles_ocean, len(pset), mean_coast, mean_coast_all, mean_ocean, mean_ocean_all]


"""
##### FUNCTION CONSTANT PUSH #####
"""
timestep_global = 5
def test_constantpush(push):
    ### START
    print "Start test_constantpush() for push={}".format(push)

    dt = timedelta(minutes=timestep_global)
    savefile = "timestep-{}".format(push)

    ### ADJUST FIELD
    fset = copy.copy(fieldset)
    coasts = stco.importCoasts("coasts-globcurrent.npz")
    stco.addCoasts(fset, coasts, constant=push)

    ### PARTICLES
    pset = ParticleSet(fieldset=fset, pclass=stc.StuckCoastParticle, lon=lons_mesh, lat=lats_mesh)

    stf.removeLandParticles(particleset=pset, filename=fieldfiles)

    ### KERNELS
    kernels = AdvectionRK4 + pset.Kernel(stc.checkVelocity) + pset.Kernel(stg.periodicBC) + pset.Kernel(stco.returnFromCoast_A) + pset.Kernel(stc.particleOnCoast)

    ### RUN
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

    ### RESULTS
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

    ### EXPORT DATA
    sta.exportParticleData(fset, pset, savefile=savefile+"-data")

    print "Stop test_constantpush for {} with {}, {}".format(push, stuckparticles, coastparticles)

    return [push, stuckparticles, coastparticles, oceanparticles, halfparticles_coast, halfparticles_ocean, len(pset), mean_coast, mean_coast_all, mean_ocean, mean_ocean_all]


"""
##### RUN TIMESTEP #####
"""

timesteps_a = np.arange(1, 30.1, 1)
timesteps_b = np.arange(5, 60.1, 5)/60.
timesteps_c = np.arange(30, 90.1, 5)
timesteps_0 = np.round(np.sort(np.unique(np.concatenate((timesteps_a, timesteps_b, timesteps_c)))), 3)
print np.size(timesteps_0)

pool = Pool(num)
results_0 = pool.map(test_timestep, timesteps_0)
results_0 = np.array(list(results_0))
np.savez_compressed("global-timestep-results", data=results_0)


"""
##### RUN PUSH #####
"""

constants_ = [0] + np.logspace(1, -5, 25).tolist() + [0]

timestep_global = 5
pool = Pool(num)
results_05 = pool.map(test_constantpush, constants_)
results_05 = np.array(list(results_05))
np.savez_compressed("global-constant-results-05", data=results_05)

timestep_global = 10
pool = Pool(num)
results_10 = pool.map(test_constantpush, constants_)
results_10 = np.array(list(results_10))
np.savez_compressed("global-constant-results-10", data=results_10)

timestep_global = 15
pool = Pool(num)
results_15 = pool.map(test_constantpush, constants_)
results_15 = np.array(list(results_15))
np.savez_compressed("global-constant-results-15", data=results_15)




print " ===== END ===== "
