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


def GlobalStuckParticlesAdvection(simulation, particles, coasts, savename, fieldlocation="GlobCurrent/", domain=None):
    """ Function for advecting particles using given parameters """
    [time_total, time_step] = simulation
    [p_lons, p_lats] = particles
    [abs, factor, constant] = coasts
    if domain is not None:
        [d_lons, d_lats] = domain
        indices = stf.getIndicesGlobCurrent(d_lons, d_lats)
    elif domain is None:
        indices = {}


    ## Create save folder
    try:
        os.makedirs(savename)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    ## Import field data
    # NOTE: Fieldset data location is hardcoded now
    fset = stf.getFieldsetGlobCurrent(filelocation=fieldlocation, time_extrapolation=True, indices=indices)

    ## Find coasts
    coasts = stf.createCoastVelocities(fieldset=fset, factor=factor, abs=abs, constant=constant)
    stf.exportCoastVelocities(coasts, savename+"coasts")

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
        output_file=pset.ParticleFile(name=savename+"trajectory", outputdt=timedelta(hours=6)),
        recovery={ErrorCode.ErrorOutOfBounds: stg.deleteParticle}
    )


    ## Export data
    return sta.exportParticleData(fieldset=fset, particleset=pset, velocities=True, savefile=savename+"particle_data")


def GlobalStuckParticlesAnalysis(parameters, savename, particledata=None, coast=False, plot=False, histogram=False, show=False, locations=False, velocities=False, trajectory=False, scatter=False, title=""):
    """"""
    [time_stuck, time_moving] = parameters

    ## Read data
    if particledata is None:
        try:
            data = sta.importParticleData(savename+"particle_data.npz")
        except:
            print "GlobalStuckParticlesAnalysis(): '{}particle_data.npz' not found".format(savename)
            return
    elif isinstance(particledata, list):
        data = particledata
    else:
        print "GlobalStuckParticlesAnalysis(): cannot find data"

    ## Process data
    subdata_all = sta.extractStuckParticles(data=data, time_stuck=0., time_moving=0., level=0, text=True)

    if coast:
        coastfield = stf.importCoastVelocities(savename+"coasts.npz")
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
        stpl.plotLocations(subdata=subdata, initial=False, show=show, savefile=savename+"final_locations", coastfields=coastfield, title=title)

    ## histogram
    if histogram:
        stpl.plotHistogram(subdata=subdata, width=5, show=show, savefile=savename+"histogram", title=title)

    ## scatter
    if scatter:
        stpl.scatterStuckMoving(subdata=subdata, show=show, savefile=savename+"scatter", title=title)

    return len(subdata)


def main():
    simulation = [75, 10] #days, minutes
    particles = [np.linspace(-175, 175, num=25), np.linspace(-75, 75, num=15)]
    coasts = [True, -1, 0]
    savename = "GlobalStuckParticles_test/"

    data = GlobalStuckParticlesAdvection(simulation, particles, coasts, savename)

    GlobalStuckParticlesAnalysis([0, 0], savename, particledata=None, coast=True, plot=True, histogram=True, show=False, locations=True, velocities=False, scatter=True)


if __name__ == "__main__":
    main()
