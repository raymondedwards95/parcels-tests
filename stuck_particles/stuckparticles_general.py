""" Some general functions used for (stuck) particles

function: haversine(x, method=None)
    np.
function: globalDistance(p1, p2)
function: absolute(x, bool=True)
kernelfunction: periodicBC(particle, fieldset, time, dt)
kernelfunction: deleteParticle(particle, fieldset, time, dt)
"""
import numpy as np


def haversine(x, method=None):
    """ default: use sin, otherwise use cos """
    if method == "cos":
        return (1. - np.cos(x)) / 2.
    else:
        return np.power(np.sin(x / 2.), 2.)


def globalDistance(p1, p2):
    """ calculate distance between two points on the earth, using haversine.
    p1 = [lon1, lat1]
    p2 = [lon2, lat2]
    * in degrees

    output = distance in m

    see https://en.wikipedia.org/wiki/Haversine_formula
    """
    [lon1, lat1] = np.radians(p1)
    [lon2, lat2] = np.radians(p2)

    r = 6371000.

    h = haversine(lat2 - lat1) + np.cos(lat1) * np.cos(lat2) * haversine(lon2 - lon1)
    return 2. * r * np.arcsin(np.power(h, 1/2.))


def absolute(x, bool=True):
    """ Returns absoulte value of the first argument if the second argument
    is true.
    """
    if bool is False:
        return x
    else:
        return np.abs(x)


def periodicBC(particle, fieldset, time, dt):
    """ Kernel used for periodic boundary conditions in east-west direction
    For correct working, also use
    `fieldset.add_periodic_halo(zonal=True)`
    """
    # from tutorials
    if particle.lon < fieldset.halo_west:
        particle.lon += fieldset.halo_east - fieldset.halo_west
    elif particle.lon > fieldset.halo_east:
        particle.lon -= fieldset.halo_east - fieldset.halo_west


def deleteParticle(particle, fieldset, time, dt):
    """ Kernel used when ErrorOutOfBounds occurs. """
    particle.delete()
    print "ErrorOutOfBounds --> Delete Particle {} at ({}, {})".format(particle.id, particle.lon, particle.lat)


def DeleteParticle(particle, fieldset, time, dt):
    """ For compatibility """
    return deleteParticle(particle, fieldset, time, dt)
