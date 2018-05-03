""" Functions and kernels for particles near coasts
function: findCoasts(filelocation=None, fieldset=None, times=[0], indices={})
    np.
    stf.getFieldsetGlobCurrent
    stf.checkCoast1d
function: exportCoasts(coastfields, filename)
function: importCoasts(filename)
function: addCoasts(fieldset, coastfields, constant=0.01)
    parcels.Field as Field
kernelfunction: returnFromCoast_A(particle, fieldset, time, dt)
    math.cos
    math.pi
kernelfunction: returnFromCoast_B(particle, fieldset, time, dt)
kernelfunction: returnFromCoast_C(particle, fieldset, time, dt)
    math.cos
    math.pi
kernelfunction: returnFromCoast_D(particle, fieldset, time, dt)
kernelfunction: returnFromCoast(particle, fieldset, time, dt)
"""
import stuckparticles_field as stf
import stuckparticles_plot as stpl

from parcels import Field

import numpy as np
import math


def findCoasts(filelocation=None, fieldset=None, times=[0], indices={}):
    """ Returns 4 2-d arrays which contain information about locations of
    coasts.
    Example: first array (N) is zero everywhere except for points that have
    ocean at the point north of it.
    """
    if filelocation is None and fieldset is None:
        print "findCoasts(): no fields found, returning"
        return
    elif filelocation is not None:
        fieldset = stf.getFieldsetGlobCurrent(filelocation, indices=indices, full_load=True)

    shape_U = fieldset.U.data.shape
    shape_V = fieldset.V.data.shape

    if not times:
        times = np.arange(shape_U[0])

    # needs clean-up:
    f_n = np.zeros(tuple((len(times),)) + shape_V[1:], dtype=np.bool)
    f_s = np.zeros(tuple((len(times),)) + shape_V[1:], dtype=np.bool)
    f_e = np.zeros(tuple((len(times),)) + shape_U[1:], dtype=np.bool)
    f_w = np.zeros(tuple((len(times),)) + shape_U[1:], dtype=np.bool)

    grid_n = grid_s = fieldset.V.grid
    grid_e = grid_w = fieldset.U.grid

    for t in times:
        # replace [0] with time

        # Start with N/S
        for i in range(f_n.shape[-1]):
            for j in range(1, f_n.shape[-2]-1):
                check = stf.checkCoast1d(fieldset, i, j, direction="y")
                if check[0]:
                    f_s[t, j, i] = True
                elif check[1]:
                    f_n[t, j, i] = True

        # Do the same for E/W
        for i in range(f_e.shape[-1]):
            for j in range(1, f_e.shape[-2]-1):
                check = stf.checkCoast1d(fieldset, i, j, direction="x")
                if check[0]:
                    f_w[t, j, i] = True
                elif check[1]:
                    f_e[t, j, i] = True

    f_n = np.mean(f_n, axis=0)
    f_s = np.mean(f_s, axis=0)
    f_e = np.mean(f_e, axis=0)
    f_w = np.mean(f_w, axis=0)

    north = ["F_N", f_n, grid_n]
    south = ["F_S", f_s, grid_s]
    east = ["F_E", f_e, grid_e]
    west = ["F_W", f_w, grid_w]

    return [north, south, east, west]


def exportCoasts(coastfields, filename):
    """ Save coastfields from findCoasts() to a file"""
    [north, south, east, west] = coastfields

    if not isinstance(filename, str):
        print "exportCoasts(): 'filename' is not a string"
        return

    np.savez_compressed(filename, north=north, south=south, east=east, west=west)
    print "exportCoasts(): fields saved as '{}'".format(filename)


def importCoasts(filename):
    """ Import coastfields exported by exportCoasts() """
    data = np.load(filename)
    return [data["north"], data["south"], data["east"], data["west"]]


def addCoasts(fieldset, coastfields, constant=0.01, show=False):
    """ Add coastfields from findCoasts() or importCoasts() and a constant to fieldset """
    fieldset.add_constant("f_constant", constant)

    for [name, data, grid] in coastfields:
        new_field = Field(name, data, lon=grid.lon, lat=grid.lat, interp_method='nearest', allow_time_extrapolation=True, transpose=False)
        fieldset.add_field(new_field)
        if show:
            print "addCoasts(): '{}' added to fieldset".format(name)


def returnFromCoast_A(particle, fieldset, time, dt):
    """ Kernel for pushing particles back from coast to ocean
    Assuming constant as 'velocity' in m/s
    Converting to lon/lat by using a simple conversion
    """
    lon, lat, depth = particle.lon, particle.lat, particle.depth

    particle.lon += dt * fieldset.f_constant / (111111 * math.cos(lat * math.pi / 180.)) * (fieldset.F_E[time, lon, lat, depth] - fieldset.F_W[time, lon, lat, depth])
    particle.lat += dt * fieldset.f_constant / 111111 * (fieldset.F_N[time, lon, lat, depth] - fieldset.F_S[time, lon, lat, depth])


def returnFromCoast_B(particle, fieldset, time, dt):
    """ Kernel for pushing particles back from coast to ocean
    Assuming constant as 'velocity' in deg/s
    """
    lon, lat, depth = particle.lon, particle.lat, particle.depth

    particle.lon += dt * fieldset.f_constant * (fieldset.F_E[time, lon, lat, depth] - fieldset.F_W[time, lon, lat, depth])
    particle.lat += dt * fieldset.f_constant * (fieldset.F_N[time, lon, lat, depth] - fieldset.F_S[time, lon, lat, depth])


def returnFromCoast_C(particle, fieldset, time, dt):
    """ Kernel for pushing particles back from coast to ocean
    Assuming constant as 'displacement' in m
    Converting to lon/lat by using a simple conversion
    """
    lon, lat, depth = particle.lon, particle.lat, particle.depth

    particle.lon += fieldset.f_constant / (111111 * math.cos(lat * math.pi / 180.)) * (fieldset.F_E[time, lon, lat, depth] - fieldset.F_W[time, lon, lat, depth])
    particle.lat += fieldset.f_constant / 111111 * (fieldset.F_N[time, lon, lat, depth] - fieldset.F_S[time, lon, lat, depth])


def returnFromCoast_D(particle, fieldset, time, dt):
    """ Kernel for pushing particles back from coast to ocean
    Assuming constant as 'displacement' in deg
    """
    lon, lat, depth = particle.lon, particle.lat, particle.depth

    particle.lon += fieldset.f_constant * (fieldset.F_E[time, lon, lat, depth] - fieldset.F_W[time, lon, lat, depth])
    particle.lat += fieldset.f_constant * (fieldset.F_N[time, lon, lat, depth] - fieldset.F_S[time, lon, lat, depth])


def returnFromCoast(particle, fieldset, time, dt):
    """ Kernel for pushing particles back from coast to ocean
    Assuming constant as 'displacement' in deg
    Used for compatibility
    This is the same as returnFromCoast_D
    """
    lon, lat, depth = particle.lon, particle.lat, particle.depth

    particle.lon += fieldset.f_constant * (fieldset.F_E[time, lon, lat, depth] - fieldset.F_W[time, lon, lat, depth])
    particle.lat += fieldset.f_constant * (fieldset.F_N[time, lon, lat, depth] - fieldset.F_S[time, lon, lat, depth])


def main():
    import stuckparticles_field as stf
    import stuckparticles_plot as stpl

    # load field using 'full_load=False'
    fset = stf.getFieldsetGlobCurrent(filelocation="GlobCurrent/")

    # find coasts (also try export and import)
    coasts = findCoasts(filelocation="GlobCurrent/")
    exportCoasts(coasts, filename="coasts-globcurrent")
    coasts = importCoasts(filename="coasts-globcurrent.npz")

    # add coasts to fieldset
    addCoasts(fset, coasts)

    # print (new) fields in fieldset
    num = len(fset.fields)
    print "fields in fieldset", fset
    for i in range(num):
        print " ", fset.fields[i].name

    stpl.plotCoasts(coasts)

if __name__ == "__main__":
    main()
