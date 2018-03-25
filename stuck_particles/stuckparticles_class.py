""" Classes and kernels for observing particles """
from parcels import FieldSet, ParticleSet, Variable, JITParticle
from parcels import AdvectionRK4, plotTrajectoriesFile, ErrorCode, ScipyParticle

import stuckparticles as st

import numpy as np
import math
import matplotlib.pyplot as plt

from datetime import timedelta, datetime
from operator import attrgetter


def periodicBC(particle, fieldset, time, dt):
    # from tutorials
    if particle.lon < fieldset.halo_west:
        particle.lon += fieldset.halo_east - fieldset.halo_west
    elif particle.lon > fieldset.halo_east:
        particle.lon -= fieldset.halo_east - fieldset.halo_west


def DeleteParticle(particle, fieldset, time, dt):
    particle.delete()
    print "ErrorOutOfBounds --> Delete Particle {} at ({}, {})".format(particle.id, particle.lon, particle.lat)


class stuckParticle(JITParticle):
    """ add more information to particles to check if a particle is stuck

    new variables:
    * time_stuck - time since last large movement/velocity

    * last_distance - last traveled distance
    * last_velocity - last velocity

    * init_lon - initial longitude (does not change)
    * init_lat - initial latitude (does not change)
    * init_time - initial time (does not change)

    * prev_lon - previous longitude
    * prev_lat - previous latitude
    * prev_time - previous time
    """
    time_stuck = Variable('time_stuck', initial=0., dtype=np.float32)
    time_moving = Variable('time_moving', initial=0., dtype=np.float32)

    last_distance = Variable('last_distance', dtype=np.float32, to_write=False, initial=0.)
    last_velocity = Variable('last_velocity', dtype=np.float32, to_write=False, initial=0.)

    init_lon = Variable('init_lon', dtype=np.float32, to_write=True, initial=attrgetter('lon'))
    init_lat = Variable('init_lat', dtype=np.float32, to_write=True, initial=attrgetter('lat'))
    init_time = Variable('init_time', dtype=np.float32, to_write=True, initial=attrgetter('time'))

    prev_lon = Variable('prev_lon', dtype=np.float32, to_write=False, initial=attrgetter('lon'))
    prev_lat = Variable('prev_lat', dtype=np.float32, to_write=False, initial=attrgetter('lat'))
    prev_time = Variable('prev_time', dtype=np.float32, to_write=False, initial=attrgetter('time'))


def checkVelocity(particle, fieldset, time, dt):
    """ Function to calculate distance and velocity and to use them to check if a particle is stuck
    For use as kernel.

    sleep_time in s
    vel_threshold in m/s

    for calculating distance we use haversine
    """
    sleep_time = 0.#5.*24*60*60 # to ensure particles move at least a few days, can also be done with time_threshold
    vel_threshold = math.pow(10., -3.)
    time_threshold = 5.*24*60*60

    r = 6371000. # radius of Earth

    if time > sleep_time:
        # prev_lon_r = math.radians(particle.prev_lon)
        # prev_lat_r = math.radians(particle.prev_lat)
        # lon_r = math.radians(particle.lon)
        # lat_r = math.radians(particle.lat)
        prev_lon_r = particle.prev_lon * math.pi / 180.
        prev_lat_r = particle.prev_lat * math.pi / 180.
        lon_r = particle.lon * math.pi / 180.
        lat_r = particle.lat * math.pi / 180.


        # note: haversine(x) = math.pow(math.sin(x / 2.), 2.)
        h = math.pow(math.sin((lat_r - prev_lat_r) / 2.), 2.) + math.cos(prev_lat_r) * math.cos(lat_r) * math.pow(math.sin((lon_r - prev_lon_r) / 2.), 2.)

        particle.last_distance = 2. * r * math.asin(math.pow(h, 1/2.))
        # delta_time = particle.time - particle.prev_time
        delta_time = dt
        particle.last_velocity = particle.last_distance / delta_time

        ## NEW
        # if particle.last_velocity < vel_threshold and particle.time_moving > time_threshold:
        #     # i.e. particle velocity is lower than threshold
        #     # and particle has moved for a certain time
        #     particle.time_stuck += delta_time
        # elif particle.last_velocity => vel_threshold
        #     # particle is moving
        #     particle.time_moving += delta_time
        #     if particle.time_stuck > 0:
        #         particle.time_stuck = 0.
        # else:
        #     # particle velocity is low and particle has not moved for a certain time
        #     if particle.time_stuck > 0:
        #         particle.time_stuck = 0.

        ## OLD
        # if particle.last_velocity < vel_threshold and particle.time_moving > time_threshold
        if particle.last_velocity < vel_threshold:
            # i.e. particle velocity is lower than threshold
            particle.time_stuck += delta_time
        else:
            particle.time_moving += delta_time
            if particle.time_stuck > 0:
                particle.time_stuck = 0.


        particle.prev_time = particle.time
        particle.prev_lon = particle.lon
        particle.prev_lat = particle.lat


def exportParticleData(fieldset, particleset, velocities=False, savefile=None):
    """ Export data per particle:
    id, lon, lat, time, init_lon, init_lat, init_time, time_stuck, time_moving

    if velocities==True: also export
    last gridvelocities
    """
    data = []
    for particle in particleset:
        sublist = np.array([particle.id,
                            particle.lon,
                            particle.lat,
                            particle.time,
                            particle.init_lon,
                            particle.init_lat,
                            particle.init_time,
                            particle.time_stuck,
                            particle.time_moving])
        if velocities:
            gridpoints = st.getGridPoints(fieldset, [particle.time, particle.lon, particle.lat, particle.depth])
            sublist = np.append(sublist, st.getGridVelocity(fieldset, gridpoints))

        data.append(sublist)

    if savefile:
        np.savez_compressed(savefile, data)
    return data






def removeLandParticles(fieldset, particleset, show=False):
    """ Remove particles that are on land.
    Particles on land are particles with velocities 0 at
    grid points around.
    """
    list_coords = st.particleCoords(particleset)
    n = 0
    for i in range(np.shape(list_coords)[0]-1, 0, -1):
        gridpoints = st.getGridPoints(fieldset, list_coords[i])
        gridvelocity = st.getGridVelocity(fieldset, gridpoints)

        # check if all u and v are 0
        if not np.any(gridvelocity[0]) and not np.any(gridvelocity[1]):
            if show:
                print "removeLandParticles(): deleted particle", particleset[i].id
            particleset.remove(i)
            n = n + 1
    print "removeLandParticles(): deleted {} particles".format(n)

    # for particle in particleset:
    #     coords = st.particleCoords(particleset, int(particle.id))
    #     gridpoints = st.getGridPoints(fieldset, coords, radius=1)
    #     gridvelocity = st.getGridVelocity(fieldset, gridpoints)
    #
    #     # check if all u and v are 0
    #     if not np.any(gridvelocity[0]) and not np.any(gridvelocity[1]):
    #         particle.delete()
    #         print "removeLandParticles(): deleted particle", particle.id
