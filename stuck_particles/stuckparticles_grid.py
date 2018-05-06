""" Additional (grid-based) particle functions for observing particles

function: getGridPoints(fieldset, coords, radius=1)
    np.
function: getGridVelocity(fieldset, vector, method="data")
function: calculateFlux(vector)
    spg.globalDistance
function: checkFlux(fluxes)
"""
import stuckparticles_general as spg

import numpy as np


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
    lats.sort(reverse=True)
    xs.sort()
    ys.sort(reverse=True)

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
        print "getGridVelocity(): recieved lons:", lons
        print "getGridVelocity(): recieved lats:", lats
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


def calculateFlux(vector):
    """ check if fluid disappears by calculating fluxes using data from getGridVelocity()
    """
    if len(vector) == 6:
        [u, v, lons, lats, depth, time] = vector
    else:
        print "calculateFlux(): vector not in correct shape. Expected 6, got", len(vector)

    res0 = np.size(lons)
    res1 = np.size(lats)
    if res0 == res1:
        res = res0
    else:
        print "calculateFlux(): Grid is not square, i.e. len(lons) != len(lats)"
        return

    if res == 2:
        [lon_l, lon_r] = lons
        [lat_u, lat_l] = lats
    else:
        half = int(np.round(res/2))
        [lon_l, lon_r] = lons[half-1 : half+1]
        [lat_u, lat_l] = lats[half-1 : half+1]
        u = u[half-1 : half+1, half-1 : half+1]
        v = v[half-1 : half+1, half-1 : half+1]

    u = u.flatten()
    v = v.flatten()

    d_north = globalDistance([lon_l, lat_u], [lon_r, lat_u])
    d_south = globalDistance([lon_l, lat_l], [lon_r, lat_l])
    d_east = globalDistance([lon_r, lat_u], [lon_r, lat_l])
    d_west = globalDistance([lon_l, lat_u], [lon_l, lat_l])

    # using average velocities
    v_north = 1/2 * (v[0] + v[1])
    v_south = 1/2 * (v[2] + v[3])
    u_east = 1/2 * (u[1] + u[3])
    u_west = 1/2 * (u[0] + u[2])

    f_north = d_north * v_north
    f_south = d_south * v_south
    f_east = d_east * u_east
    f_west = d_west * u_west

    return [f_north, f_south, f_east, f_west]


def checkFlux(fluxes):
    """ Check if total flux is close to 0 """
    [f_north, f_south, f_east, f_west] = fluxes
    return f_north - f_south + f_east - f_west <= np.power(10., -6.)
