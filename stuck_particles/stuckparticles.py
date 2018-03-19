from parcels import FieldSet, ParticleSet, Variable, JITParticle
from parcels import AdvectionRK4, plotTrajectoriesFile, ErrorCode, ScipyParticle

import numpy as np
import math
import matplotlib.pyplot as plt

from datetime import datetime


def DeleteParticle(particle, fieldset, time, dt):
    p = particle
    p.delete()
    print "ErrorOutOfBounds --> Delete Particle {} at ({}, {})".format(p.id, p.lon, p.lat)


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

        return [times, lons, lats, depths]

    elif type(identification) == int:
        for i in range(len(particleset)):
            p = particleset[i]
            if particleset[i].id == identification:
                times = p.time
                lons = p.lon
                lats = p.lat
                depths = p.depth

                return [times, lons, lats, depths]

    elif type(identification) == list:
        for i in range(len(particleset)):
            p = particleset[i]
            if particleset[i].id in identification:
                times.append(p.time)
                lons.append(p.lon)
                lats.append(p.lat)
                depths.append(p.depth)

        return [times, lons, lats, depths]

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


def plotVelocity(vector, field="U", coords=None, show=None, savefile=None, vmax=None):
    """ plot results of getVelocity() """
    if savefile is None:
        show = True
    if show is None and savefile is not None:
        show = False

    if len(vector) == 6:
        [U, V, lons, lats, depth, time] = vector
    else:
        print "plotVelocity(): vector not in correct shape. Expected 6, got", len(vector)
        return

    if coords is not None and np.size(coords) == 4:
        [time, lon, lat, depth] = coords
        plon, plat = lon, lat
    else:
        plon = lons[len(lons)/2]
        plat = lats[len(lats)/2]

    if field == "U":
        vel = U
    elif field == "V":
        vel = V
    else:
        print "plotVelocity(): parameter 'field' is incorrect. Expected 'U' or 'V'"

    if vmax is not None and type(vmax) == float:
        vmin = -vmax
    else:
        vmin = None
        vmin = None

    if math.isnan(time):
        # at the start, time is nan,
        # so this is an ad-hoc solution
        time = "start"

    plt.figure()
    plt.contourf(lons, lats, vel, vmin=vmin, vmax=vmax)
    plt.plot(plon, plat, 'ro')
    plt.title("field {} at t={}".format(field, time))
    plt.xlabel("longitude")
    plt.ylabel("latitude")
    plt.colorbar()
    plt.grid()
    if show:
        plt.show()
    if savefile is not None:
        plt.savefig(savefile)



def plotAbsoluteVelocity(vector, coords=None, show=None, savefile=None, vmax=None):
    """ Plot results of absoluteVelocity() """
    if savefile is None:
        show = True
    if show is None and savefile is not None:
        show = False

    if len(vector) == 5:
        [vel, lons, lats, depth, time] = vector
    else:
        print "plotAbsoluteVelocity(): vector not in correct shape. Expected 5, got", len(vector)
        return

    if coords is not None and np.size(coords) == 4:
        [time, lon, lat, depth] = coords
        plon, plat = lon, lat
    else:
        plon = lons[len(lons)/2]
        plat = lats[len(lats)/2]

    if math.isnan(time):
        # at the start, time is nan,
        # so this is an ad-hoc solution
        time = "start"

    plt.figure()
    plt.contourf(lons, lats, vel, vmin=0, vmax=vmax)
    plt.plot(plon, plat, 'ro')
    plt.title("absolute velocities at t={}".format(time))
    plt.xlabel("longitude")
    plt.ylabel("latitude")
    plt.colorbar()
    plt.grid()
    if show:
        plt.show()
    if savefile is not None:
        plt.savefig(savefile)


def getGridPoints(fieldset, coords, radius=1):
    """ Return the closest grid points (lon, lat) """
    radius = int(radius)

    if np.size(coords) == 4:
        [time, lon, lat, depth] = coords
    else:
        return
        print "getGridPoints(): Arguments (coords) not in correct shape. Expected 4, got", np.size(coords)

    grid_lon = fieldset.U.grid.lon
    grid_lat = fieldset.U.grid.lat

    if lon < np.min(grid_lon) or lon > np.max(grid_lon):
        print "getGridPoints(): Longitude in not in the domain. {} is not between {} and {}.".format(lon, np.min(grid_lon), np.max(grid_lon))
        return
    if lat < np.min(grid_lat) or lat > np.max(grid_lat):
        print "getGridPoints(): Latitude in not in the domain. {} is not between {} and {}.".format(lat, np.min(grid_lat), np.max(grid_lat))
        return

    lons = []
    lats = []
    xs = []
    ys = []

    for i in range(np.size(grid_lon)-1):
        if lon < grid_lon[i]:
            for ii in range(radius):
                lons.append(grid_lon[i-ii-1])
                lons.append(grid_lon[i+ii])
                xs.append(i-ii-1)
                xs.append(i+ii)
            break

    for j in range(np.size(grid_lat)-1):
        if lat < grid_lat[j]:
            for jj in range(radius):
                lats.append(grid_lat[j-jj-1])
                lats.append(grid_lat[j+jj])
                ys.append(j-jj-1)
                ys.append(j+jj)
            break

    # could be optimized by checking if sort is needed
    lons.sort()
    lats.sort()
    xs.sort()
    ys.sort()

    return [lons, lats, xs, ys, depth, time]


def getGridVelocity(fieldset, vector, method="data"):
    """ Use results from getGridPoints() to find velocities on given points
    method="data" or method="eval"
    """
    if len(vector) == 6:
        [lons, lats, xs, ys, depth, time] = vector
    else:
        print "getGridVelocity(): vector not in correct shape. Expected 6, got", len(vector)

    res0 = np.size(lons)
    res1 = np.size(lats)
    res2 = np.size(xs)
    res3 = np.size(ys)
    if res0 == res1 == res2 == res3:
        res = res0
    else:
        print "getGridVelocity(): Grid is not square, i.e. len(lons) != len(lats)"
        return

    u = np.zeros((res, res))
    v = np.zeros((res, res))

    if method == "data":
        for i in range(res):
            for j in range(res):
                u[i][j] = fieldset.U.data[0, ys[i], xs[j]]
                v[i][j] = fieldset.V.data[0, ys[i], xs[j]]
    elif method == "eval":
        for i in range(res):
            for j in range(res):
                u[j][i] = fieldset.U.eval(time, lons[i], lats[j], depth)
                v[j][i] = fieldset.V.eval(time, lons[i], lats[j], depth)
    else:
        print "getGridVelocity(): method not found, use \"data\" or \"eval\"."
        return

    return [u, v, lons, lats, depth, time]


def main():
    pass



if __name__ == "__main__":
    main()
