""" Tools for analysing data for stuck particles """
import stuckparticles as st
import stuckparticles_class as stc

import numpy as np
import math

try:
    import matplotlib.pyplot as plt
except:
    pass


def exportParticleData(fieldset, particleset, velocities=False, savefile=None):
    """ Export data per particle:
    id, lon, lat, time, init_lon, init_lat, init_time, time_stuck, time_moving, time_simulated

    if velocities == true: also export
    grid_u, grid_v, grid_lons, grid_lats, grid_depth, grid_time
    """
    data = []

    for p in particleset:
        sublist = np.array([p.id, p.lon, p.lat, p.time, p.init_lon, p.init_lat, p.init_time, p.time_stuck, p.time_moving, p.time_simulated])

        if velocities:
            gp = st.getGridPoints(fieldset, [p.time, p.lon, p.lat, p.depth])
            gv = st.getGridVelocity(fieldset, gp)
            # sublist.append(gv)
            sublist = np.append(sublist, gv)

        data.append(sublist)

    if savefile:
        np.savez_compressed(savefile, data)
        print "Data saved in", savefile
    return data


def importParticleData(filename):
    """ Read data exported using exportParticleData() """
    return np.load(filename)["arr_0"]


def extractStuckParticles(data, time_stuck=5, time_moving=5, level=0, text=False):
    """ Get data for particles that are stuck for a number of days after moving for a number of days.

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
        print "Total number of particles: {}. {} are 'stuck', {} are 'free'.".format(len(data), stuck_particles, free_particles)
    return list


def printLocations(subdata, initial=False):
    """ Show initial and last location of particles in 'subdata' """
    if initial and subdata[0][-1] < 3:
        print "printLocations(): missing initial lon and lat."
        initial = False

    if initial:
        for p in subdata:
            print "Particle {:.0f}: Initial location ({:05.3f}, {:05.3f}) --> Last location ({:05.3f}, {:05.3f}).".format(p[0], p[5], p[6], p[1], p[2])

    else:
        for p in subdata:
            print "Particle {:.0f}: Last location ({:05.3f}, {:05.3f}).".format(p[0], p[1], p[2])



def plotLocations(subdata, title="", initial=False, show=None, savefile=None):
    """ Plot locations of particles in subdata or data
    Assuming that the first three values are (id, lon, lat).
    """
    if savefile is None:
        show = True
    if show is None and savefile is not None:
        show = False

    if initial and subdata[0][-1] < 3:
        print "plotLocations(): not enough data."
        initial = False

    colors = ["red", "blue", "green", "pink", "orange", "purple"]
    number = len(colors)

    plt.figure()
    for i in range(len(subdata)):
        m = i%number

        plt.plot(subdata[i][1], subdata[i][2], "o", markersize=3, color=colors[m], label="Particle {} at ({:.1f}, {:.1f})".format(int(subdata[i][0]), subdata[i][1], subdata[i][2]))

        if initial and subdata[0][-1] > 2:
            plt.plot(subdata[i][5], subdata[i][6], "o", markersize=1, color=colors[m])
            plt.plot([subdata[i][1], subdata[i][5]], [subdata[i][2], subdata[i][6]], "--", color=colors[m])

    plt.grid()
    plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    plt.xlabel("longitude")
    plt.ylabel("latitude")
    plt.title(title)

    if show:
        plt.show()
    if savefile is not None:
        plt.savefig(savefile)


def printGridVelocity(subdata, flux=False, index=None):
    """ Show U and V on the grid points around the particles
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
                flux_numbers = st.calculateFlux(p[7:13])
                flux_check = st.checkFlux(flux_numbers)
                print "Flux on grid:"
                print flux_numbers
                print flux_check
