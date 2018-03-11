from parcels import FieldSet, ParticleSet, Variable, JITParticle
from parcels import AdvectionRK4, plotTrajectoriesFile, ErrorCode, ScipyParticle

import numpy as np
import math
import matplotlib.pyplot as plt

def DeleteParticle(particle, fieldset, time, dt):
    p = particle
    p.delete()
    print "ErrorOutOfBounds --> Delete Particle {} at ({}, {})".format(p.id, p.lon, p.lat)

def particleCoords(particleset, identification=None):
    """ returns an array with time, lon, lat and depth of
    particles in the particleset with ids in the list
    identification.
    """
    coords = []
    if identification == None:
        for i in range(len(particleset)):
            p = particleset[i]
            coords.append([p.time, p.lon, p.lat, p.depth])
        return np.array(coords)
    if type(identification) == list:
        for i in range(len(particleset)):
            p = particleset[i]
            if particleset[i].id in identification:
                coords.append([p.time, p.lon, p.lat, p.depth])
        return np.array(coords)
    return np.array(coords)

def getVelocity(fieldset, coords, radius=0, step=1):
    """ return U, V for lons and lats close to given coords with
    form (time, lon, lat, depth)
    radius and step is used to also
    get U, V at (lon,lat) +-step with step <= radius
    result is an array of size 1+2*radius/step
    """
    if np.size(coords) == 4:
        lon, lat = coords[0][1], coords[0][2]
        time, depth = coords[0][0], coords[0][3]
    else:
        print "coords not in correct shape"
        return

    if step == 0 and radius == 0:
        print "* step and radius are 0, setting step to 1"
        step = 1
    elif step == 0 and radius != 0:
        print "* step is 0, setting step equal to radius,", radius
        step = radius
    if step > radius and radius != 0:
        print "* step is larger than radius, setting step equal to radius,", radius
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

    return [U, V, lons, lats]


def absoluteVelocity(vector):
    """ Calculate absolute velocities from list [U, V, lon, lat]: (U^2 + v^2)^(1/2) """
    if len(vector) == 4:
        U, V = vector[0], vector[1]
        lon, lat = vector[2], vector[3]
    else:
        print "vector not in correct shape"

    res0 = np.shape(U)[0]
    res1 = np.shape(U)[1]
    res2 = np.shape(V)[0]
    res3 = np.shape(V)[1]
    if res0 == res1 == res2 == res3:
        res = res0
    else:
        print "U and V are not square or do not have the same shapes."
        return

    vel = np.ones((res, res))

    for i in range(res):
        for j in range(res):
            vel[i, j] = np.sqrt(np.power(U[i, j], 2) + np.power(V[i, j], 2))

    return [vel, lon, lat]


def plotAbsoluteVelocity(vector, savefile=None, vmax=3):
    """ Plot results of absoluteVelocity() """
    if len(vector) == 3:
        vel = vector[0]
        lon = vector[1]
        lat = vector[2]
        plon = lon[len(lon)/2]
        plat = lat[len(lat)/2]
    else:
        print "vector not in correct shape"
        return
    plt.figure()
    plt.contourf(lon, lat, vel, vmin=0, vmax=vmax)
    plt.plot(plon, plat)
    plt.xlabel("longitude")
    plt.ylabel("latitude")
    plt.colorbar()
    if savefile is None:
        plt.show()
    else:
        plt.savefig(savefile)


def getGridPoints(fieldset, coords, radius=1):
    """ Return the closest grid points (lon, lat) """
    if np.size(coords) == 4:
        lon, lat = coords[0][1], coords[0][2]
        time, depth = coords[0][0], coords[0][3]
    else:
        print "coords not in correct shape"
        return

    grid_lon = fieldset.U.grid.lon
    grid_lat = fieldset.U.grid.lat

    if lon < np.min(grid_lon) or lon > np.max(grid_lon):
        print "Longitude in not in the domain"
        return
    if lat < np.min(grid_lat) or lat > np.max(grid_lat):
        print "Latitude in not in the domain"
        return

    for i in range(np.size(grid_lon)-1):
        if lon < grid_lon[i]:
            lon_0 = grid_lon[i-1]
            lon_1 = grid_lon[i]
            x_0 = i-1
            x_1 = i
            break
    for j in range(np.size(grid_lat)-1):
        if lat < grid_lat[j]:
            lat_0 = grid_lat[j-1]
            lat_1 = grid_lat[j]
            y_0 = j-1
            y_1 = j
            break

    lons = [lon_0, lon_1]
    lats = [lat_0, lat_1]
    xs = [x_0, x_1]
    ys = [y_0, y_1]

    return [lons, lats, xs, ys]


def getGridVelocity(fieldset, vector):
    """ Use results from getGridPoints() to find velocities on given points """
    [lons, lats, xs, ys] = vector

    u = np.zeros((2, 2))
    v = np.zeros((2, 2))

    for i in range(2):
        for j in range(2):
            u[i][j] = fieldset.U.data[0, ys[i], xs[j]]
            v[i][j] = fieldset.V.data[0, ys[i], xs[j]]

    return [u, v, lons, lats]


def main():
    pass


if __name__ == "__main__":
    main()
