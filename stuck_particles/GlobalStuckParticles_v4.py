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


pool = pool(2)

timesteps_a = np.arange(1, 30.1, 1)
timesteps_b = np.arange(5, 60.1, 5)/60.
timesteps_c = np.arange(30, 90.1, 5)
timesteps_0 = np.flip(np.round(np.sort(np.unique(np.concatenate((timesteps_a, timesteps_b, timesteps_c)))), 3), -1)
print np.size(timesteps_0)

results_0 = pool.map(timestepTest, timesteps_0)
results_0 = np.array(list(results_0))
np.savez_compressed("global-timestep-results", data=results_0)


print " ===== END ===== "
