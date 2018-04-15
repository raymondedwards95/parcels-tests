""" Field manipulation functions for stuck particles

function: checkField(fieldset, text=True)
    np.
function: getIndicesGlobCurrent(lons, lats)
function: getFieldsetGlobCurrent(filelocation, indices={}, time_extrapolation=False)
    parcels.FieldSet
function: checkCoast1d(fieldset, x, y, direction=None, time=0)
function: createCoastVelocities(fieldset, factor=True, abs=True, constant=0)
    stg.absolute
function: addGlobCurrentCoast(fieldset, coastfields)
function: exportCoastVelocities(field_coast_U, field_coast_V, filename)
function: importCoastVelocities(filename)
function: removeLandParticles(fieldset, particleset, show=False)
    stp.particleCoords
    stgr.getGridPoints
    stgr.getGridVelocity

function: main()
    stpl.showCoast
"""
from parcels import FieldSet

import numpy as np

import stuckparticles_particles as stp
import stuckparticles_general as stg
import stuckparticles_grid as stgr


def checkField(fieldset, text=True):
    """ Check if lons, lats are the same for different fields in fieldset """
    if text:
        print "\nFieldset contains the following fields:"
        for i in range(len(fieldset.fields)):
            print fieldset.fields[i].name

    ulon = fieldset.U.grid.lon
    ulat = fieldset.U.grid.lat
    udep = fieldset.U.grid.depth
    vlon = fieldset.V.grid.lon
    vlat = fieldset.V.grid.lat
    vdep = fieldset.V.grid.depth

    if text:
        if np.all(ulon == vlon):
            print "longitudes are the same for U and V"
        else:
            print "longitudes are not the same for U and V. Note that not all functions will work as intended."
        if np.all(ulat == vlat):
            print "latitudes are the same for U and V"
        else:
            print "latitudes are not the same for U and V. Note that not all functions will work as intended."
        if np.all(udep == vdep):
            print "depths are the same for U and V"
        else:
            print "depths are not the same for U and V. Note that not all functions will work as intended."

    return np.all(ulon == vlon) and np.all(ulat == vlat) and np.all(udep == vdep)


def getIndicesGlobCurrent(lons, lats):
    """ Find indices for GlobCurrent data """
    if np.size(lons) == 1:
        lon_0, lon_1 = int(np.floor(lons-5)), int(np.ceil(lons+5))
    else:
        lon_0, lon_1 = int(np.round(np.min(lons))), int(np.round(np.max(lons)))

    if np.size(lats) == 1:
        lat_0, lat_1 = int(np.floor(lats-5)), int(np.ceil(lats+5))
    else:
        lat_0, lat_1 = int(np.round(np.min(lats))), int(np.round(np.max(lats)))

    lon_range = range((lon_0-5+180)*4-1, (lon_1+5+180)*4+1)
    lat_range = range((lat_0-5+80)*4-1, (lat_1+5+80)*4+1)

    indices = {"lon": lon_range,
               "lat": lat_range}

    print "getIndicesGlobCurrent(): Success! Indices created."
    return indices


def getFieldsetGlobCurrent(filelocation, indices={}, time_extrapolation=False):
    """ Create fieldset from GlobCurrent data """
    filenames = {"U": filelocation+"*.nc",
                 "V": filelocation+"*.nc"}
    var = {"U": "eastward_eulerian_current_velocity",
           "V": "northward_eulerian_current_velocity"}
    dim = {"lat": "lat",
           "lon": "lon",
           "time": "time"}
    fieldset = FieldSet.from_netcdf(filenames, variables=var, dimensions=dim, indices=indices, allow_time_extrapolation=time_extrapolation)

    print "getFieldsetGlobCurrent(): Success! Fieldset imported."
    return fieldset


def checkCoast1d(fieldset, x, y, direction=None, time=0):
    """ Check if the given point-indices (x, y) is an coast point.
    Method is by checking for [0, 0, U] or for [U, 0, 0], with center point as given point.
    Use direction ("x" or "y") to find coasts perpendicular to the given direction. If None, then check for both "x" and "y".
    For "x" or "y" return two bools:
    Use index time to choose a different field-time.
    First return is for [ocean, coast, land] and the second is for [land, coast, ocean].
    """
    if direction == None:
        coast_x = checkCoast1d(fieldset, x, y, direction="x")
        coast_y = checkCoast1d(fieldset, x, y, direction="y")
        return coast_x, coast_y

    elif direction == "x":
        vector_U = fieldset.U.data[time, y, x-1:x+2]
        vector_U_trim = np.trim_zeros(vector_U)
        if len(vector_U_trim) == 1:
            # Checks if vector contains 2 zeros and one non-zero
            # and if the non_zero is at the begin or end
            if vector_U_trim == vector_U[0]:
                return [True, False]
            elif vector_U_trim == vector_U[-1]:
                return [False, True]
            else:
                return [False, False]
        else:
            return [False, False]

    elif direction == "y":
        vector_V = fieldset.V.data[time, y-1:y+2, x]
        vector_V_trim = np.trim_zeros(vector_V)
        if len(vector_V_trim) == 1:
            # Checks if vector contains 2 zeros and one non-zero
            # and if the non_zero is at the begin or end
            if vector_V_trim == vector_V[0]:
                return [True, False]
            elif vector_V_trim == vector_V[-1]:
                return [False, True]
            else:
                return [False, False]
        else:
            return [False, False]

    else:
        print "checkCoast1d(): direction is not None, 'x' or 'y'."
        return False


def createCoastVelocities(fieldset, factor=True, abs=True, constant=0):
    """ Search for gridpoints that represent coasts, i.e. groups of
    points with velocities [U, 0, 0] or [0, 0, U] (or V).
    Then calculate the new center value as +- abs(U * factor) +- constant.

    Returns two fields (U and V) with 'imaginary'-velocities at the coast.

    Function works for fieldsets with fields U and V with
    two spatial coordinates (lon, lat) and one temporal coordinate (time).
    """
    vel_U = fieldset.U.data
    vel_V = fieldset.V.data
    dims_U = vel_U.shape
    dims_V = vel_V.shape
    lons_U = fieldset.U.grid.lon
    lons_V = fieldset.V.grid.lon
    lats_U = fieldset.U.grid.lat
    lats_V = fieldset.V.grid.lat

    if factor == 0 and constant == 0:
        factor = True

    if len(dims_U) == 3 and len(dims_V) == 3:
        if dims_U[-1] == dims_V[-1]:
            nx = dims_U[-1]
        else:
            print "createCoastVelocities(): shapes of U.data and V.data are not the same."
            return

        if dims_U[-2] == dims_V[-2]:
            ny = dims_U[-2]
        else:
            print "createCoastVelocities(): shapes of U.data and V.data are not the same."
            return
    else:
        print "createCoastVelocities(): dimensions of U.data and V.data are not 3, but {} and {}.".format(len(dims_U), len(dims_V))
        return

    if np.any(lons_U != lons_V) or np.any(lats_U != lats_V):
        print "createCoastVelocities(): grids of U and V are not the same."
        return
    else:
        lons = lons_U
        lats = lats_U

    field_coast_U = np.zeros([ny, nx])
    field_coast_V = np.zeros([ny, nx])

    for x in range(0, nx):
        # print "createCoastVelocities(): current index {} of {}.".format(x, nx)
        # Do for all longitudes, using modulo (%)
        for y in range(0+1, ny-1):
        # Do for all latitudes, but not top or bottom
            # Check in x direction
            check_x = checkCoast1d(fieldset, x, y, direction="x")
            if check_x[0]:
                field_coast_U[y, x] = -constant + -stg.absolute(factor * vel_U[0, y, (x-1) % nx], abs)
            elif check_x[1]:
                field_coast_U[y, x] = constant + stg.absolute(factor * vel_U[0, y, (x+1) % nx], abs)

            # Check in y direction
            check_y = checkCoast1d(fieldset, x, y, direction="y")
            if check_y[0]:
                field_coast_V[y, x] = -constant + -stg.absolute(factor * vel_V[0, y-1, x], abs)
            elif check_y[1]:
                field_coast_V[y, x] = constant + stg.absolute(factor * vel_V[0, y+1, x], abs)

    if factor == True:
        field_coast_U = field_coast_U.astype(np.bool)
        field_coast_V = field_coast_V.astype(np.bool)

    return [field_coast_U, field_coast_V, lons, lats]


def addGlobCurrentCoast(fieldset, coastfields):
    """ Returns a new fieldset that results from
    fieldset + coastfields
    """
    [field_coast_U, field_coast_V, lons, lats] = coastfields
    new_fieldset = fieldset

    if np.shape(fieldset.U.data)[0] == np.shape(fieldset.V.data)[0]:
        nt = np.shape(fieldset.U.data)[0]
    else:
        print "addGlobCurrentCoast(): fieldset.U.data and fieldset.V.data do not have the same shape."
        return
    if field_coast_U.dtype == bool and field_coast_U.dtype == bool:
        pass
    else:
        for t in range(nt):
            new_fieldset.U.data[t] += field_coast_U
            new_fieldset.V.data[t] += field_coast_V

    return new_fieldset


def exportCoastVelocities(coastfields, filename):
    """ Export fields as numpy-array-files """
    [field_coast_U, field_coast_V, lons, lats] = coastfields

    if not isinstance(filename, str):
        print "exportCoastVelocities(): filename is not a string"
        return
    np.savez_compressed(filename, coast_U=field_coast_U, coast_V=field_coast_V, lons=lons, lats=lats)
    print "exportCoastVelocities(): Fields saved as " + filename

def importCoastVelocities(filename):
    data = np.load(filename)
    return [data["coast_U"], data["coast_V"], data["lons"], data["lats"]]


def removeLandParticles(fieldset, particleset, show=False):
    """ Remove particles that are on land.
    Particles on land are particles with velocities 0 at
    grid points around.
    """
    print "removeLandParticles(): there are {} particles in the ParticleSet".format(len(particleset))
    list_coords = stp.particleCoords(particleset)
    n = 0
    for i in range(np.shape(list_coords)[0]-1, 0, -1):
        gridpoints = stgr.getGridPoints(fieldset, list_coords[i])
        gridvelocity = stgr.getGridVelocity(fieldset, gridpoints)

        # check if all u and v are 0
        if not np.any(gridvelocity[0]) and not np.any(gridvelocity[1]):
            if show:
                print "removeLandParticles(): deleted particle", particleset[i].id
            particleset.remove(i)
            n = n + 1
    print "removeLandParticles(): deleted {} particles that were on land".format(n)
    print "removeLandParticles(): there are now {} particles in the ParticleSet".format(len(particleset))


def main():
    import stuckparticles_plot as stpl

    fset = getFieldsetGlobCurrent("GlobCurrent/", time_extrapolation=True)
    coast_fields = createCoastVelocities(fset)
    exportCoastVelocities(coast_fields[0], coast_fields[1], "GlobCurrentCoast")

    new_fset = addGlobCurrentCoast(fset, coast_fields)

    stpl.showCoast(coast_fields)


if __name__ == "__main__":
    main()
