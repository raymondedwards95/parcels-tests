""" Particle functions for observing particles

function: particleCoords(particleset, identification=None)
    np.
function: getVelocity(fieldset, coords, radius=0, step=1)
function: absoluteVelocity(vector)

"""
import numpy as np


def particleCoords(particleset, identification=None):
    """ returns an array with time, lon, lat and depth of
    particles in the particleset with ids in
    identification (int, list or None).
    """
    times = []
    lons = []
    lats = []
    depths = []

    if identification == None:
        for i in range(len(particleset)):
            p = particleset[i]
            times.append(p.time)
            lons.append(p.lon)
            lats.append(p.lat)
            depths.append(p.depth)

        return np.transpose([times, lons, lats, depths]).tolist()

    elif isinstance(identification, (int, np.integer)):
        for i in range(len(particleset)):
            p = particleset[i]
            if particleset[i].id == identification:
                times = p.time
                lons = p.lon
                lats = p.lat
                depths = p.depth

                return [times, lons, lats, depths]

    elif isinstance(identification, (list)):
        for i in range(len(particleset)):
            p = particleset[i]
            if particleset[i].id in identification:
                times.append(p.time)
                lons.append(p.lon)
                lats.append(p.lat)
                depths.append(p.depth)

        return np.transpose([times, lons, lats, depths]).tolist()

    else:
        print "particleCoords(): 'indentification' (", identification, ") is incorrect"

    print "particleCoords(): Something went wrong, returning data from first particle"
    p = particleset[0]
    return [p.time, p.lon, p.lat, p.depth]


def getVelocity(fieldset, coords, radius=0, step=1):
    """ return U, V for lons and lats close to given coords of one particle, with
    form (time, lon, lat, depth)
    radius and step is used to also
    get U, V at (lon,lat) +-step with step <= radius
    result is an array of size 1+2*radius/step
    """
    if np.size(coords) == 4:
        [time, lon, lat, depth] = coords
    else:
        print "getVelocity(): coords not in correct shape. Expected 4, got", np.size(coords)
        return

    if step == 0 and radius == 0:
        print "getVelocity(): step and radius are 0, setting step to 1"
        step = 1
    elif step == 0 and radius != 0:
        print "getVelocity(): step is 0, setting step equal to radius,", radius
        step = radius
    if step > radius and radius != 0:
        print "getVelocity(): step is larger than radius, setting step equal to radius,", radius
        step = radius

    res = int(1+2*np.round(radius/step))

    lons = np.linspace(lon-radius, lon+radius, res)
    lats = np.linspace(lat-radius, lat+radius, res)

    U = np.ones((res, res))
    V = np.ones((res, res))

    for i in range(res):
        for j in range(res):
            U[j][i] = fieldset.U.eval(time, lons[i], lats[j], depth)
            V[j][i] = fieldset.V.eval(time, lons[i], lats[j], depth)

    return [U, V, lons, lats, depth, time]



def absoluteVelocity(vector):
    """ Calculate absolute velocities from list [U, V, lons, lats, time]: (U^2 + v^2)^(1/2) """
    if len(vector) == 6:
        [U, V, lons, lats, depth, time] = vector
    else:
        print "absoluteVelocity(): vector not in correct shape. Expected 6, got", len(vector)

    res0 = np.shape(U)[0]
    res1 = np.shape(U)[1]
    res2 = np.shape(V)[0]
    res3 = np.shape(V)[1]
    if res0 == res1 == res2 == res3:
        res = res0
    else:
        print "absoluteVelocity(): U and V are not square or do not have the same shapes. Shapes of U and V are", np.shape(U), np.shape(V)
        return

    vel = np.ones((res, res))

    for i in range(res):
        for j in range(res):
            vel[i, j] = np.sqrt(np.power(U[i, j], 2) + np.power(V[i, j], 2))

    return [vel, lons, lats, depth, time]
