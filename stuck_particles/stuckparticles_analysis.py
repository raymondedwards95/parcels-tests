""" Tools for analysing data for stuck particles

function: exportParticleData(fieldset, particleset, velocities=False, savefile=None)
    np.
    stgr.getGridPoints
    stgr.getGridVelocity
function: importParticleData(filename)
# function: extractStuckParticles(data, time_stuck=5, time_moving=5, level=0, text=False)
function: def rearrangeData(data, level=0)
function: printLocations(subdata, initial=False)
# function: printGridVelocity(subdata, flux=False, index=None)
#     stgr.calculateFlux
#     stgr.checkFlux
function: filterStuckParticles(subdata, current_time_stuck=None, current_time_moving=None, total_time_stuck=None, total_time_moving=None)
function: filterCoastParticles(subdata, current_time_coast=None, current_time_ocean=None, total_time_coast=None, total_time_ocean=None)
"""
import stuckparticles_grid as stgr

import numpy as np
import math


def exportParticleData(fieldset, particleset, velocities=False, savefile=None, StuckParticle=True, CoastParticle=True):
    """ Export data per particle:

    * StuckParticle: id, lon, lat, time, init_lon, init_lat, init_time, time_stuck, time_moving, total_time_stuck, total_time_moving, time_simulated
    * CoastParticle: id, lon, lat, time, total_time_coast, total_time_ocean, current_time_coast, current_time_ocean, number_on_coast, number_in_ocean, time_simulated
    * stuckCoastParticle: id, lon, lat, time, init_lon, init_lat, init_time, time_stuck, time_moving, total_time_stuck, total_time_moving, total_time_coast, total_time_ocean, current_time_coast, current_time_ocean, number_on_coast, number_in_ocean, time_simulated

    if velocities == true: also export
    grid_u, grid_v, grid_lons, grid_lats, grid_depth, grid_time
    """
    if len(particleset) == 0:
        print "exportParticleData(): no particles in 'particleset', returning"
        return

    # determine class of particles
    p = particleset[0]
    if StuckParticle:
        try:
            if p.time_stuck >= 0.:
                c_StuckParticle = True
                print "exportParticleData(): 'particleset' contains StuckParticle-data"
            else:
                c_StuckParticle = False
                print "exportParticleData(): 'particleset' does not contain StuckParticle-data"
        except:
            c_StuckParticle = False
            print "exportParticleData(): 'particleset' does not contain StuckParticle-data!"

    if CoastParticle:
        try:
            if p.total_time_coast >= 0.:
                c_CoastParticle = True
                print "exportParticleData(): 'particleset' contains CoastParticle-data"
            else:
                c_CoastParticle = False
                print "exportParticleData(): 'particleset' does not contain CoastParticle-data"
        except:
            c_CoastParticle = False
            print "exportParticleData(): 'particleset' does not contain CoastParticle-data!"


    data = []

    for p in particleset:
        if c_StuckParticle and not c_CoastParticle:
            sublist = np.array([p.id, p.lon, p.lat, p.time, p.init_lon, p.init_lat, p.init_time, p.time_stuck, p.time_moving, p.total_time_stuck, p.total_time_moving, p.time_simulated])
        elif not c_StuckParticle and c_CoastParticle:
            sublist = np.array([p.id, p.lon, p.lat, p.time, p.total_time_coast, p.total_time_ocean, p.current_time_coast, p.current_time_ocean, p.number_on_coast, p.number_in_ocean, p.time_simulated])
        elif c_StuckParticle and c_CoastParticle:
            sublist = np.array([p.id, p.lon, p.lat, p.time, p.init_lon, p.init_lat, p.init_time, p.time_stuck, p.time_moving, p.total_time_stuck, p.total_time_moving, p.total_time_coast, p.total_time_ocean, p.current_time_coast, p.current_time_ocean, p.number_on_coast, p.number_in_ocean, p.time_simulated])
        else:
            sublist = np.array([p.id, p.lon, p.lat, p.time])

        if velocities:
            gp = stgr.getGridPoints(fieldset, [p.time, p.lon, p.lat, p.depth])
            gv = stgr.getGridVelocity(fieldset, gp)
            sublist = np.append(sublist, gv)

        data.append(sublist)

    if savefile:
        np.savez_compressed(savefile, particledata=data)
        print "Data saved in", savefile
    return data


def importParticleData(filename):
    """ Read data exported using exportParticleData() """
    return np.load(filename)["particledata"]


def extractStuckParticles(data, time_stuck=5, time_moving=5, level=0, text=False):
    """ OLD - Get data for particles that are stuck for a number of days after moving for a number of days.
    For use with specific particle class: StuckParticle.

    Number of elements depends on 'level':
    0 - everything
    1 - simple (id, lon, lat)
    2 - basic (id, lon, lat, time_stuck, time_moving)
    3 - history (id, lon, lat, time_stuck, time_moving, init_lon, init_lat)
    4 - gridinformation (id, lon, lat, time_stuck, time_moving, init_lon, init_lat, grid_u, grid_v, grid_lons, grid_lats, grid_depth, grid_time)
    """
    stuck_particles = 0
    free_particles = 0

    if level == 0:
        ## For compability, set level to max if level == 0
        level = 4

    list = []

    for p in data:
        if p[7] >= time_stuck*24*60*60 and p[8] >= time_moving*24*60*60:
            stuck_particles += 1
            if level == 0:
                sublist = p.tolist()
                sublist.append(level)
                list.append(sublist)

            elif level > 0:
                sublist = p[0:3].tolist()

                if level > 1:
                    sublist.append(p[7])
                    sublist.append(p[8])

                if level > 2:
                    sublist.append(p[4])
                    sublist.append(p[5])

                if level > 3:
                    sublist.append(p[10])
                    sublist.append(p[11])
                    sublist.append(p[12])
                    sublist.append(p[13])
                    sublist.append(p[14])
                    sublist.append(p[15])

                sublist.append(level)
                list.append(sublist)

    if text:
        print "extractStuckParticles(): Total number of particles: {}. Using given arguments: {} are 'stuck' for at least {} days, {} are 'free'.".format(len(data), stuck_particles, time_stuck, free_particles)
    return list


def rearrangeData(data, level=0, show=False):
    """ Reduce the amount of data, depending on given 'level'.
    Note: all a-lists will come before all b-lists
    Note 2: a-lists are only possible if class of particles was StuckParticle
            b-lists are only possible if class was CoastParticle

    0 - 'everything'
    1 - id, lon, lat <--> [3], [4], [5]

    2a - total_time_stuck, total_time_moving <--> [6], [7]
    2b - total_time_coast, total_time_ocean <--> [-1], [-2]

    3a - current_time_stuck, current_time_moving <--> [8], [9]
    3b - current_time_coast, current_time_ocean <--> [-3], [-4]

    4b - number_on_coast, number_in_ocean <--> [-5], [-6]

    5a - init_lon, init_lat, init_time <--> [10], [11], [12]

    For level = 0:
    return [level, bool:StuckParticle, bool:CoastParticle, id, lon, lat, total_time_stuck, total_time_moving, current_time_stuck, current_time_moving, init_lon, init_lat, init_time, number_in_ocean, number_on_coast, current_time_ocean, current_time_coast, total_time_ocean, total_time_coast] for each particle
    """
    if level <= 0 or level > 5:
        level = 5
        if show:
            print "rearrangeData(): Setting level to highest possible:", level

    if len(data[0]) == 12:
        c_StuckParticle = True
        c_CoastParticle = False
    elif len(data[0]) == 11:
        c_StuckParticle = False
        c_CoastParticle = True
    elif len(data[0]) == 18:
        c_StuckParticle = True
        c_CoastParticle = True
    elif len(data[0]) == 4:
        c_StuckParticle = False
        c_CoastParticle = False
    else:
        print "rearrangeData(): length of data is incorrect, expected 4, 11, 12 or 18 instead of", len(data[0])
        return

    list = []

    for p in data:
        sublist = [level, c_StuckParticle, c_CoastParticle]

        if level >= 1:
            sublist.append(p[0])
            sublist.append(p[1])
            sublist.append(p[2])
        else:
            print "rearrangeData(): ERROR, level (={}) is not larger than 0. Returning".format(level)
            return

        if c_StuckParticle is True:
            if level >= 2:
                sublist.append(p[9])
                sublist.append(p[10])
            if level >= 3:
                sublist.append(p[7])
                sublist.append(p[8])
            if level >= 4:
                pass
            if level >= 5:
                sublist.append(p[4])
                sublist.append(p[5])
                sublist.append(p[6])

        if c_CoastParticle is True:
            if level >= 5:
                pass
            if level >= 4:
                sublist.append(p[-3])
                sublist.append(p[-2])
            if level >= 3:
                sublist.append(p[-4])
                sublist.append(p[-5])
            if level >= 2:
                sublist.append(p[-6])
                sublist.append(p[-7])

        list.append(sublist)

    if show:
        if c_StuckParticle:
            print "rearrangeData(): data contains 'StuckParticle'-data up to and including 'level' {}".format(level)
        if c_CoastParticle:
            print "rearrangeData(): data contains 'StuckParticle'-data up to and including 'level' {}".format(level)

    return list

    ##### * StuckParticle: id, lon, lat, time, init_lon, init_lat, init_time, time_stuck, time_moving, total_time_stuck, total_time_moving, time_simulated
    ##### * CoastParticle: id, lon, lat, time, total_time_coast, total_time_ocean, current_time_coast, current_time_ocean, number_on_coast, number_in_ocean, time_simulated
    ##### * stuckCoastParticle: id, lon, lat, time, init_lon, init_lat, init_time, time_stuck, time_moving, total_time_coast, total_time_ocean, current_time_coast, current_time_ocean, number_on_coast, number_in_ocean, time_simulated


def printLocations(subdata, initial=False, show=True):
    """ Show and return initial and last location of particles in 'subdata' """
    if initial and subdata[0][0] < 5:
        print "printLocations(): missing initial lon and lat."
        initial = False

    locs = []

    if initial:
        for p in subdata:
            if show:
                print "Particle {:.0f}: Initial location ({:05.3f}, {:05.3f}) --> Last location ({:05.3f}, {:05.3f}).".format(p[3], p[10], p[11], p[4], p[5])
            locs.append([p[3], p[4], p[5], p[10], p[11]])

    else:
        for p in subdata:
            if show:
                print "Particle {:.0f}: Last location ({:05.3f}, {:05.3f}).".format(p[3], p[4], p[5])
            locs.append([p[3], p[4], p[5]])

    if show:
        print

    return locs


def printGridVelocity(subdata, flux=False, index=None):
    """ OLD - Show U and V on the grid points around the particles
    if flux==True: calculate flux in the gridcell
    """
    if subdata[0][-1] < 4:
        print "printLocations(): missing grid information."
        return

    if index == None:
        index = np.arange(0, len(subdata), dtype=np.int32).tolist()
    elif isinstance(index, (int, np.integer)):
        index = [index]

    for i in range(len(subdata)):
        if i in index:
            p = subdata[i]
            print "\nParticle {:.0f} at ({:05.3f}, {:05.3f}):".format(p[0], p[1], p[2])
            print "U on grid:"
            print p[7]
            print "V on grid:"
            print p[8]

            if flux:
                flux_numbers = stgr.calculateFlux(p[7:13])
                flux_check = stgr.checkFlux(flux_numbers)
                print "Flux on grid:"
                print flux_numbers
                print flux_check
    print ""


def filterStuckParticles(subdata, current_time_stuck=None, current_time_moving=None, total_time_stuck=None, total_time_moving=None):
    """ Function to find particles with (more/different) specific parameters
    using particles in 'subdata', created by rearrangeData().
    Filters are in days.

    Function returns two lists, one with 'wanted' particles and one with 'unwanted' particles
    """
    DAYS = 24.*60.*60.

    list_a = []
    list_b = []

    if subdata[0][0] < 2:
        print "filterStuckParticles(): can't filter using 'total_time_stuck', 'total_time_moving', 'current_time_stuck' or 'current_time_moving', not enough information in 'subdata'. Returning an empty and a full list/"
        for p in subdata:
            list_b.append(p)
        return (list_a, list_b)
    if subdata[0][0] < 3:
        print "filterStuckParticles(): can't filter using 'current_time_stuck' or 'current_time_moving', not enough information. Continuing for 'total_time_stuck' and 'total_time_moving'."
        current_time_moving = current_time_stuck = None

    if subdata[0][1] is False:
        print "filterStuckParticles(): can't filter using 'total_time_stuck', 'total_time_moving', 'current_time_stuck' or 'current_time_moving', class of particle is not 'StuckParticle'. Returning an empty and a full list"
        for p in subdata:
            list_b.append(p)
        return (list_a, list_b)

    # if current_time_stuck is not None or current_time_moving is not None:
    #     if current_time_stuck is None:
    #         current_time_stuck = 0
    #     elif current_time_moving is None:
    #         current_time_moving = 0
    # if total_time_stuck is not None or total_time_moving is not None:
    #     if total_time_stuck is None:
    #         total_time_stuck = 0
    #     elif total_time_moving is None:
    #         total_time_moving = 0
    if current_time_stuck is None:
        current_time_stuck = 0.
    if current_time_moving is None:
        current_time_moving = 0.
    if total_time_stuck is None:
        total_time_stuck = 0.
    if total_time_moving is None:
        total_time_moving = 0.

    if subdata[0][0] >= 3:
        for p in subdata:
            if p[3] > current_time_stuck*DAYS and p[4] > current_time_moving*DAYS and p[6] > total_time_stuck*DAYS and p[7] > total_time_moving*DAYS:
                list_a.append(p)
            else:
                list_b.append(p)
    else:
        for p in subdata:
            if p[6] > total_time_stuck*DAYS and p[7] > total_time_moving*DAYS:
                list_a.append(p)
            else:
                list_b.append(p)


    print "filterStuckParticles(): New number of particles: {} and {}. ".format(len(list_a), len(list_b))

    return (list_a, list_b)


def filterCoastParticles(subdata, current_time_coast=None, current_time_ocean=None, total_time_coast=None, total_time_ocean=None):
    """ Function to find particles with (more/different) specific parameters
    using particles in 'subdata', created by rearrangeData().
    Filters are in days.

    Function returns two lists, one with 'wanted' particles and one with 'unwanted' particles
    """
    DAYS = 24.*60.*60.

    list_a = []
    list_b = []

    if subdata[0][0] < 2:
        print "filterCoastParticles(): can't filter using 'total_time_coast', 'total_time_ocean', 'current_time_coast' or 'current_time_ocean', not enough information in 'subdata'. Returning an empty and a full list/"
        for p in subdata:
            list_b.append(p)
        return (list_a, list_b)
    if subdata[0][0] < 3:
        print "filterCoastParticles(): can't filter using 'current_time_coast' or 'current_time_ocean', not enough information. Continuing for 'total_time_coast' and 'total_time_ocean'."
        current_time_moving = current_time_stuck = None

    if subdata[0][2] is False:
        print "filterCoastParticles(): can't filter using 'total_time_coast', 'total_time_ocean', 'current_time_coast' or 'current_time_ocean', class of particle is not 'CoastParticle'. Returning an empty and a full list"
        for p in subdata:
            list_b.append(p)
        return (list_a, list_b)

    if current_time_coast is None:
        current_time_coast = 0.
    if current_time_ocean is None:
        current_time_ocean = 0.
    if total_time_coast is None:
        total_time_coast = 0.
    if total_time_ocean is None:
        total_time_ocean = 0.

    if subdata[0][0] >= 3:
        for p in subdata:
            if p[3] > current_time_coast*DAYS and p[4] > current_time_ocean*DAYS and p[6] > total_time_coast*DAYS and p[7] > total_time_ocean*DAYS:
                list_a.append(p)
            else:
                list_b.append(p)
    else:
        for p in subdata:
            if p[6] > total_time_coast*DAYS and p[7] > total_time_ocean*DAYS:
                list_a.append(p)
            else:
                list_b.append(p)


    print "filterCoastParticles(): New number of particles: {} and {}. ".format(len(list_a), len(list_b))

    return (list_a, list_b)


### [level, bool:StuckParticle, bool:CoastParticle, id, lon, lat, total_time_stuck, total_time_moving, current_time_stuck, current_time_moving, init_lon, init_lat, init_time, number_in_ocean, number_on_coast, current_time_ocean, current_time_coast, total_time_ocean, total_time_coast]
