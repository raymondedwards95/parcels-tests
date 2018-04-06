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


def GlobalStuckParticlesAdvection(simulation, particles, coasts, filename):
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
    stf.exportCoastVelocities(coasts[0], coasts[1], filename+"coasts")

    ## Add coasts to fieldset
    fset = stf.addGlobCurrentCoast(fieldset=fset, coastfields=coasts)


    ## Create particles
    p_lon_mesh, p_lat_mesh = np.meshgrid(p_lons, p_lats)
    pset = ParticleSet(fieldset=fset, pclass=stc.stuckParticle, lon=p_lon_mesh, lat=p_lat_mesh)

    ## Remove landparticles
    stf.removeLandParticles(fieldset=fset, particleset=pset)


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


def GlobalStuckParticlesAnalysis(parameters, filename, particledata=None, coast=False, plot=False, histogram=False, show=False, locations=False, velocities=False):
    """"""
    [time_stuck, time_moving] = parameters

    ## Read data
    if particledata is None:
        try:
            data = sta.importParticleData(filename+"particle_data.npz")
        except:
            print "GlobalStuckParticlesAnalysis(): '{}particle_data' not found".format(filename)
            return
    elif isinstance(particledata, list):
        data = particledata
    else:
        print "GlobalStuckParticlesAnalysis(): cannot find data"

    ## Process data
    subdata_all = sta.extractStuckParticles(data=data, time_stuck=0., time_moving=0., level=0, text=True)

    if coast:
        coastfield = stf.importCoastVelocities(filename+"coasts.npz")
    else:
        coastfield = None

    ## Filter data
    subdata = sta.filterParticles(subdata_all, time_stuck=time_stuck, time_moving=time_moving)

    ## Print locations/velocities
    if velocities:
        sta.printGridVelocity(subdata=subdata, flux=True)
    elif locations:
        sta.printLocations(subdata=subdata, initial=True)

    ## Plots
    if plot:
        stpl.plotLocations(subdata=subdata, title="Final locations of stuck particles", initial=False, show=show, savefile=filename+"final_locations", coastfield=coastfield)

    ## histogram
    if histogram:
        stpl.plotHistogram(subdata=subdata, width=1, show=show, savefile=filename+"histogram")


def main():
    simulation = [50, 10] #days, minutes
    particles = [np.linspace(-175, 175, num=11), np.linspace(-75, 75, num=9)]
    coasts = [True, -1, 0]
    filename = "GlobalStuckParticles_test/"

    data = GlobalStuckParticlesAdvection(simulation, particles, coasts, filename)

    GlobalStuckParticlesAnalysis([0, 0], filename, particledata=None, coast=True, plot=True, histogram=True, show=False, locations=True, velocities=False)


if __name__ == "__main__":
    main()
