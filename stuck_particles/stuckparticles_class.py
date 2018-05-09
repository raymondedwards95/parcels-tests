""" Classes and kernels for observing particles

class: StuckParticle
    parcels.JITParticle as JITParticle
    parcels.ScipyParticle as ScipyParticle
    parcels.Variable as Variable
    np.
    operator.attrgetter as attrgetter
kernelfunction: checkVelocity(particle, fieldset, time, dt)
    math.
class: CoastParticle
kernelfunction: particleOnCoast(particle, fieldset, time, dt)
class: StuckCoastParticle
    StuckParticle
    CoastParticle
class: stuckParticle
    StuckParticle
"""
from parcels import JITParticle, ScipyParticle, Variable

import numpy as np
import math

from operator import attrgetter


class StuckParticle(JITParticle):
    """ add more information to particles to check if a particle is stuck

    new variables:
    * time_stuck - time since last large movement/velocity
    * time_moving - time since last small movement/velocity
    * time_simulated - total time since start (should be equal to time of particle)

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
    time_simulated = Variable('time_simulated', initial=0., dtype=np.float32)

    last_distance = Variable('last_distance', dtype=np.float32, to_write=False, initial=0.)
    last_velocity = Variable('last_velocity', dtype=np.float32, to_write=False, initial=0.)

    init_lon = Variable('init_lon', dtype=np.float32, to_write=True, initial=attrgetter('lon'))
    init_lat = Variable('init_lat', dtype=np.float32, to_write=True, initial=attrgetter('lat'))
    init_time = Variable('init_time', dtype=np.float32, to_write=True, initial=attrgetter('time'))

    prev_lon = Variable('prev_lon', dtype=np.float32, to_write=False, initial=attrgetter('lon'))
    prev_lat = Variable('prev_lat', dtype=np.float32, to_write=False, initial=attrgetter('lat'))
    prev_time = Variable('prev_time', dtype=np.float32, to_write=False, initial=attrgetter('time'))

    # particle_class = Variable('particle_class', initial="StuckParticle", dtype=np.str, to_write=False)


def checkVelocity(particle, fieldset, time, dt):
    """ Function to calculate distance and velocity and to use them to check if a particle is stuck
    For use as kernel for the class StuckParticle.

    sleep_time in s
    vel_threshold in m/s

    for calculating distance we use haversine
    """
    vel_threshold = math.pow(10., -3.)

    r = 6371000. # radius of Earth

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

    if particle.last_velocity < vel_threshold:
        # i.e. particle velocity is lower than threshold
        particle.time_stuck += delta_time
    else:
        particle.time_moving += delta_time
        if particle.time_stuck > 0:
            particle.time_stuck = 0.
            particle.time_moving = 0.

    particle.prev_time = particle.time
    particle.prev_lon = particle.lon
    particle.prev_lat = particle.lat

    particle.time_simulated += dt


class CoastParticle(JITParticle):
    """ Particle class to track if a particle gets on the coast

    new variables:
    * total_time_coast - total time on coast
    * current_time_coast - current time on coast

    * total_time_ocean - total time in ocean
    * current_time_ocean - current time in ocean

    * number_on_coast - number of times on coast (ocean -> coast)
    * number_in_ocean - number of times released from coast (coast -> ocean)

    * time_simulated - total time since start (should be equal to time of particle)
    """
    total_time_coast = Variable('total_time_coast', initial=0., dtype=np.float32, to_write=True)
    current_time_coast = Variable('current_time_coast', initial=0., dtype=np.float32, to_write=True)

    total_time_ocean = Variable('total_time_ocean', initial=0., dtype=np.float32, to_write=True)
    current_time_ocean = Variable('current_time_ocean', initial=0., dtype=np.float32, to_write=True)

    number_on_coast = Variable('number_on_coast', initial=0, dtype=np.int32, to_write=True)
    number_in_ocean = Variable('number_in_ocean', initial=0, dtype=np.int32, to_write=True)
    initial_location = Variable('initial_location', initial=0, dtype=np.int32, to_write=True)

    time_simulated = Variable('time_simulated', initial=0., dtype=np.float32, to_write=True)

    # particle_class = Variable('particle_class', initial="CoastParticle", dtype=np.str, to_write=False)


def particleOnCoast(particle, fieldset, time, dt):
    """ Function to use as kernel for the class CoastParticle
    and coastfield F_N, F_S, F_E, F_W """

    # check if particle on coast:
    if (fieldset.F_N[time, particle.lon, particle.lat, particle.depth] or
            fieldset.F_S[time, particle.lon, particle.lat, particle.depth] or
            fieldset.F_E[time, particle.lon, particle.lat, particle.depth] or
            fieldset.F_W[time, particle.lon, particle.lat, particle.depth]):

        if particle.number_on_coast == 0 and particle.number_in_ocean == 0:
            particle.initial_location = 1

        if particle.current_time_coast == 0.:
            particle.number_on_coast += 1

        particle.total_time_coast += dt
        particle.current_time_coast += dt

        if particle.current_time_ocean > 0:
            particle.current_time_ocean = 0.

    else:
        if particle.number_on_coast == 0 and particle.number_in_ocean == 0:
            particle.initial_location = 0

        if particle.current_time_ocean == 0.:
            particle.number_in_ocean += 1

        particle.total_time_ocean += dt
        particle.current_time_ocean += dt

        if particle.current_time_coast > 0:
            particle.current_time_coast = 0.

    particle.time_simulated += dt


class StuckCoastParticle(StuckParticle, CoastParticle):
    """ Combine both StuckParticle and CoastParticle classes """

    time_simulated = Variable('time_simulated', initial=0., dtype=np.float32)

    # particle_class = Variable('particle_class', initial="StuckCoastParticle", dtype=np.str, to_write=False)


class stuckParticle(StuckParticle):
    """ For compatibility """

    time_simulated = Variable('time_simulated', initial=0., dtype=np.float32)
