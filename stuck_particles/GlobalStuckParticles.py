""" Global stuck particles """
from parcels import ParticleSet, AdvectionRK4, ErrorCode

import os
import errno
import numpy as np

from datetime import timedelta

import stuckparticles_analysis as sta
import stuckparticles_class as stc
import stuckparticles_field as stf
import stuckparticles_general as stg
import stuckparticles_grid as stgr
import stuckparticles_particles as stp
import stuckparticles_plot as stpl


def testGlobalStuckParticlesAdvect(simulation, particles, coasts, filename):
    """ Function for advecting particles using given parameters """
    [time_total, time_step] = simulation
    [p_lons, p_lats] = particles
    [abs, factor, constant] = coasts

    ## Create save folder
    try:
        os.makedirs(filename)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    ## Import field data
    # NOTE: Fieldset data location is hardcoded now
    fset = stf.getFieldsetGlobCurrent(filelocation="GlobCurrent/", time_extrapolation=True)

    ## Find coasts
    coasts = stf.createCoastVelocities(fieldset=fset, factor=factor, abs=abs, constant=constant)

    ## Add coasts to fieldset
    fset = stf.addGlobCurrentCoast(fieldset=fset, coastfields=coasts)


    ## Create particles
    p_lon_mesh, p_lat_mesh = np.meshgrid(p_lons, p_lats)
    pset = ParticleSet(fieldset=fset, pclass=stc.stuckParticle, lon=p_lon_mesh, lat=p_lat_mesh)

    ## Remove landparticles
    stf.removeLandParticles(fieldset=fset, particleset=pset)
    # should print number of particles


    ## Define kernels
    ## Advection
    kernels = AdvectionRK4

    ## Velocity
    k_velocity = pset.Kernel(stc.checkVelocity)
    kernels += k_velocity

    ## Periodic boundary conditions
    fset.add_periodic_halo(zonal=True)
    k_boundary = pset.Kernel(stg.periodicBC)
    kernels += k_boundary


    ## Advect
    pset.execute(
        kernels,
        runtime=timedelta(days=time_total),
        dt=timedelta(minutes=time_step),
        output_file=pset.ParticleFile(name=filename+"trajectory", outputdt=timedelta(hours=3)),
        recovery={ErrorCode.ErrorOutOfBounds: stg.deleteParticle}
    )


    ## Export data
    return sta.exportParticleData(fieldset=fset, particleset=pset, velocities=True, savefile=filename+"particle_data")
