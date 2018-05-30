""" Some general functions used for (stuck) particles

function: haversine(x, method=None)
    np.
function: globalDistance(p1, p2)
function: absolute(x, bool=True)
kernelfunction: periodicBC(particle, fieldset, time, dt)
kernelfunction: deleteParticle(particle, fieldset, time, dt)
"""
import numpy as np
import bisect


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


def sort_col(array, index=[0]):
    """ Sort array by using columns
    https://stackoverflow.com/a/2828371
    """
    if len(array.shape) != 2:
        print "sort_col(): array is not 2d, returning unsorted array"
        return array

    num = array.shape[1]
    shape = 'f' + (num-1)*',f'

    type = array.dtype

    ord = []
    for i in index:
        ord.append('f'+str(i))

    return np.sort(array.view(shape), order=ord, axis=0).view(type)


def nearest_index(a, x):
    i = bisect.bisect(a, x)
    if i == np.size(a): i -= 1
    return i


def main():
    x = np.array([[2,2,2], [1,1,1], [3,3,3]])
    y = sort_col(x, [1])
    print x
    print y
    print

    x = np.random.randint(0, 10, (5,5)).astype(np.float32)
    y = sort_col(x, [1])
    print x
    print y

    u = np.arange(0, 100, 2)
    print nearest_index(u, 3)


if __name__ == '__main__':
    main()
