""" Functions to make plots for stuck particles

plotfunction: plotVelocity(vector, field="U", coords=None, show=None, savefile=None, vmax=None)
    np.
    math.
    plt.
plotfunction: plotAbsoluteVelocity(vector, coords=None, show=None, savefile=None, vmax=None)
plotfunction: showCoast(fields, type=np.bool, show=None, savefile=None, field="all")
plotfunction: plotLocations(subdata, title="", initial=False, show=None, savefile=None, coastfields=None, coasttype=np.bool)
plotfunction: plotHistogram(subdata, width=1, show=None, savefile=None, title="")
plotfunction: scatterStuckMoving(subdata, show=None, savefile=None, title="")
plotfunction: plotCoast(coasts, show=None, savefile=None)
plotfunction: plotTrajectories(filename, ocean_particles=True, coast_particles=True, coasts=None)
    Field
    FieldSet
"""
import stuckparticles_analysis as sta

from netCDF4 import Dataset

from parcels import Field, FieldSet

import numpy as np
import math

try:
    import matplotlib.pyplot as plt
except:
    print "ERROR: matplotlib not found"


def plotVelocity(vector, field="vector", coords=None, show=None, savefile=None, vmax=None):
    """ plot results of stp.getVelocity() """
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

    mesh_x, mesh_y = np.meshgrid(lons, lats)

    if field == "U":
        vel = U
    elif field == "V":
        vel = V
    elif field == "vector":
        pass
    elif field == "absolute":
        vel = np.sqrt(np.power(U, 2) + np.power(V, 2))
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

    if field == "U" or field == "V" or field == "absolute":
        plt.contourf(lons, lats, vel, vmin=vmin, vmax=vmax)
    elif field == "vector":
        color = np.hypot(U, V)
        plt.quiver(lon, lat, U, V, color)

    plt.plot(plon, plat, 'o', markersize=5, color="red")
    plt.plot(mesh_x.flatten(), mesh_y.flatten(), 'o', markersize=3, color="black")
    # plt.title("field {} at t={}".format(field, time))
    plt.xlabel("longitude")
    plt.ylabel("latitude")
    plt.colorbar()
    plt.grid()
    if savefile is not None:
        plt.savefig(savefile)
        print "plotVelocity(): plot saved as '{}'".format(savefile)
    if show:
        print "plotVelocity(): showing plot"
        plt.show()


def plotAbsoluteVelocity(vector, coords=None, show=None, savefile=None, vmax=None):
    """ Plot results of stp.absoluteVelocity() """
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

    mesh_x, mesh_y = np.meshgrid(lons, lats)

    if math.isnan(time):
        # at the start, time is nan,
        # so this is an ad-hoc solution
        time = "start"

    plt.figure()
    plt.contourf(lons, lats, vel, vmin=0, vmax=vmax)
    plt.plot(plon, plat, 'o', markersize=5, color="red")
    plt.plot(mesh_x.flatten(), mesh_y.flatten(), 'o', markersize=3, color="black")
    plt.title("absolute velocities at t={}".format(time))
    plt.xlabel("longitude")
    plt.ylabel("latitude")
    plt.colorbar()
    plt.grid()
    if savefile is not None:
        plt.savefig(savefile)
        print "plotAbsoluteVelocity(): plot saved as '{}'".format(savefile)
    if show:
        print "plotAbsoluteVelocity(): showing plot"
        plt.show()


def showCoast(coastfields, coasttype=np.bool, show=None, savefile=None, field="all"):
    """ Plot values in fields in a contour plot """
    if len(coastfields) == 4:
        [field_coast_U, field_coast_V, lons, lats] = coastfields
    else:
        print "showCoast(): 'coastfields' contains not enough information"
        return

    if savefile is None:
        show = True
    if show is None and savefile is not None:
        show = False

    if field == "U":
        field = field_coast_U
    elif field == "V":
        field = field_coast_V
    else:
        field = np.sqrt(np.square(field_coast_U) + np.square(field_coast_V))

    field = field.astype(coasttype)

    plt.figure()
    plt.contourf(lons, lats, field)
    # plt.xlim([np.min(lons), np.max(lons)])
    # plt.ylim([np.max(lats), np.min(lats)])
    plt.title("Locations of coasts")
    plt.xlabel("longitude")
    plt.ylabel("latitude")
    plt.colorbar()
    plt.grid()
    if savefile is not None:
        plt.savefig(savefile)
        print "showCoast(): plot saved as '{}'".format(savefile)
    if show:
        print "showCoast(): showing plot"
        plt.show()


def plotLocations(subdata, title="", initial=False, show=None, savefile=None, coastfields=None, legend=False, coasttype=np.bool, filter_stuck=False, filter_coast=False):
    """ Plot locations of particles in subdata or data from stuckparticles_analysis.
    Assuming that the first three values are (id, lon, lat).
    """
    if savefile is None:
        show = True
    if show is None and savefile is not None:
        show = False

    if initial and subdata[0][0] < 5:
        print "plotLocations(): not enough data. Continuing without initial positions"
        initial = False

    if coastfields is not None:
        if len(coastfields) == 4:
            [field_coast_U, field_coast_V, lons, lats] = coastfields
        else:
            print "plotLocations(): 'coastfields' contains not enough information, skipping coasts"
            coastfields = None

        field = np.sqrt(np.square(field_coast_U) + np.square(field_coast_V))
        field = field.astype(coasttype)


    filter = None

    if filter_stuck is True and subdata[0][0] >= 3 and subdata[0][1] is True:
        filter = "filter_stuck"
        print "plotLocations(): base color of particles on state of 'time_stuck'."

    if filter_coast is True and subdata[0][0] >= 3 and subdata[0][2] is True:
        # override filter_stuck
        filter = "filter_coast"
        print "plotLocations(): base color of particles on state of 'time_coast'."


    colors = ["red", "blue", "green", "pink", "orange", "purple", "black"]
    number = len(colors)

    plt.figure()
    for i in range(len(subdata)):
        if filter == "filter_stuck":
            # now stuck: red
            # stuck before: orange
            # never stuck: green
            # always stuck: black
            if subdata[i][8] > 0. and subdata[i][7] > 0.:
                m = 0
            elif subdata[i][9] > 0. and subdata[i][6] > 0.:
                m = 4
            elif subdata[i][6] <= 0.:
                m = 2
            elif subdata[i][7] <= 0.:
                m = 6
            else:
                print "plotLocations(): 'filter_stuck' is currently not working for particle", subdata[i][3]
                m = 6

        elif filter == "filter_coast":
            # now on coast: red
            # on coast before: orange
            # never on coast: green
            # always on coast: black
            if subdata[i][-3] > 0. and subdata[i][-2] > 0.:
                m = 0
            elif subdata[i][-4] > 0. and subdata[i][-1] > 0.:
                m = 4
            elif subdata[i][-1] <= 0.:
                m = 2
            elif subdata[i][-2] <= 0.:
                m = 6
            else:
                print "plotLocations(): 'filter_coast' is currently not working for particle", subdata[i][3]
                m = 6
        else:
            m = i%number # for random colors


        plt.plot(subdata[i][4], subdata[i][5], "o", markersize=4, color=colors[m], label="Particle {} at ({:.1f}, {:.1f})".format(int(subdata[i][3]), subdata[i][4], subdata[i][5]))

        if initial and subdata[0][1] >= 5:
            plt.plot(subdata[i][10], subdata[i][11], "o", markersize=1, color=colors[m])
            plt.plot([subdata[i][4], subdata[i][10]], [subdata[i][5], subdata[i][11]], "--", color=colors[m])

    if coastfields is not None:
        plt.contourf(lons, lats, field, alpha=0.5, cmap="Greys")
        # plt.xlim([np.min(lons), np.max(lons)])
        # plt.ylim([np.max(lats), np.min(lats)])

    plt.grid()
    if legend:
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    plt.xlabel("longitude")
    plt.ylabel("latitude")
    plt.title(title)

    if savefile is not None:
        plt.savefig(savefile)
        print "plotLocations(): plot saved as '{}'".format(savefile)
    if show:
        print "plotLocations(): showing plot"
        plt.show()


def plotHistogram(subdata, width=1, show=None, savefile=None, title="", filter="coast"):
    """ Show number of stuck/coast particles in a histogram using data created in
    stuckparticles_analysis.
    Parameter filter can be "coast" or "stuck".
    """
    if savefile is None:
        show = True
    if show is None and savefile is not None:
        show = False

    if subdata[0][0] < 2:
        print "plotHistogram(): missing 'total_time' information. Exiting"
        return

    if filter == "stuck" and subdata[0][1] is False:
        print "plotHistogram(): particle class is not 'StuckParticle'. Exiting"
        return

    if filter == "coast" and subdata[0][2] is False:
        print "plotHistogram(): particle class is not 'CoastParticle'. Exiting"
        return


    list_a = np.array([])
    list_b = np.array([])

    if filter == "stuck":
        for p in subdata:
            list_a = np.append(list_a, p[6]/(24*60*60))
            list_b = np.append(list_b, p[7]/(24*60*60))
    elif filter == "coast":
        for p in subdata:
            list_a = np.append(list_a, p[-1]/(24*60*60))
            list_b = np.append(list_b, p[-2]/(24*60*60))
    else:
        print "plotHistogram(): can not use given 'filter'. Exiting"


    list_a_reduced = np.trim_zeros(sorted(list_a, reverse=True))

    # n_bins = int(round((np.max(time_stuck+time_moving) - np.min(time_stuck+time_moving)) / width))
    n_bins = np.arange(0, np.max(list_a+list_b)+width, width)

    plt.figure()
    plt.hist(list_a_reduced, n_bins, facecolor="green", alpha=0.75, edgecolor='black', linewidth=0.1)
    plt.xticks(n_bins)
    if filter == "stuck":
        plt.xlabel("days stuck")
    if filter == "coast":
        plt.xlabel("days on coast")
    plt.ylabel("number of particles")
    plt.title(title+"\n{} particles in histogram".format(len(time_stuck)))
    plt.grid()

    if savefile is not None:
        plt.savefig(savefile)
        print "plotHistogram(): plot saved as '{}'".format(savefile)
    if show:
        print "plotHistogram(): showing plot"
        plt.show()


def scatterParticleData(subdata, show=None, savefile=None, title="", filter="coast"):
    """ Create a scatter plot from time_moving/time_ocean and time_stuck/time_coast in subdata
    Plot shows if a particle was stuck/on_coast more than once.
    Parameter filter can be "coast" or "stuck".
    """
    if savefile is None:
        show = True
    if show is None and savefile is not None:
        show = False

    if subdata[0][-1] < 2:
        print "scatterStuckMoving(): missing time information."
        return

    list_a = np.array([])
    list_b = np.array([])

    if filter == "stuck":
        for p in subdata:
            list_a = np.append(list_a, p[6]/(24*60*60))
            list_b = np.append(list_b, p[7]/(24*60*60))
    elif filter == "coast":
        for p in subdata:
            list_a = np.append(list_a, p[-1]/(24*60*60))
            list_b = np.append(list_b, p[-2]/(24*60*60))
    else:
        print "scatterParticleData(): can not use given 'filter'. Exiting"

    x = np.linspace(0, np.max(list_b+list_a), 2)
    y = np.max(list_b+list_a) - x

    plt.figure()
    plt.scatter(list_a, list_b)
    plt.plot(x, y, "--", color="red")
    if filter == "stuck":
        plt.xlabel("days stuck")
        plt.ylabel("days moved")
    if filter == "coast":
        plt.xlabel("days on coast")
        plt.ylabel("days in ocean")
    plt.xlim([-0.1, max(list_b+list_a)*1.05])
    plt.ylim([-0.1, max(list_b+list_a)*1.05])
    plt.title(title+"\n{} particles in plot".format(len(subdata)))
    plt.grid()


    if savefile is not None:
        plt.savefig(savefile)
        print "scatterStuckMoving(): plot saved as '{}'".format(savefile)
    if show:
        print "scatterStuckMoving(): showing plot"
        plt.show()


def plotCoasts(coasts, show=None, savefile=None):
    if savefile is None:
        show = True
    if show is None and savefile is not None:
        show = False

    [north, south, east, west] = coasts

    # add data
    field_N = north[1].astype(np.bool)
    field_S = south[1].astype(np.bool)
    field_E = east[1].astype(np.bool)
    field_W = west[1].astype(np.bool)

    plt.figure()

    # should change to contourf
    plt.contour(north[2].lon, north[2].lat, field_N)
    plt.contour(south[2].lon, south[2].lat, field_S)
    plt.contour(east[2].lon, east[2].lat, field_E)
    plt.contour(west[2].lon, west[2].lat, field_W)

    plt.title("Locations of coasts")
    plt.xlabel("longitude")
    plt.ylabel("latitude")
    # plt.colorbar()
    plt.grid()


    if savefile is not None:
        plt.savefig(savefile)
        print "plotCoasts(): plot saved as '{}'".format(savefile)
    if show:
        print "plotCoasts(): showing plot"
        plt.show()


# def readPlotTrajectoriesFile(filename):
#     file = Dataset(filename, "r")
#     lon = np.ma.filled(file.variables["lon"], np.nan)
#     lat = np.ma.filled(file.variables["lat"], np.nan)
#     time = np.ma.filled(file.variables["time"], np.nan)
#     z = np.ma.filled(file.variables["z"], np.nan)


def plotTrajectories(filename, ocean_particles=True, coast_particles=True, field=None, coasts=None, show=None, savefile=None):
    """ Plot trajectories of particles in a trajectories-file (*.nc)
    ### to do: filters
    """
    if savefile is None:
        show = True
    if show is None and savefile is not None:
        show = False

    file = Dataset(filename, "r")

    lon = np.ma.filled(file.variables["lon"], np.nan)
    lat = np.ma.filled(file.variables["lat"], np.nan)
    time = np.ma.filled(file.variables["time"], np.nan)
    z = np.ma.filled(file.variables["z"], np.nan)

    plt.figure()
    plt.plot(np.transpose(lon), np.transpose(lat), "--")
    plt.plot(np.transpose(lon)[-1], np.transpose(lat)[-1], "o")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")

    if coasts is not None:
        if isinstance(coasts, Field):
            # coasts is a Field from a FieldSet
            plt.contourf(coasts.grid.lon, coasts.grid.lat, coasts.data[0].astype(np.bool), alpha=0.5, cmap="Greys")

        elif isinstance(coasts, (list, np.ndarray)):
            if np.shape(np.array(coasts)) == (4, 3):
                # coasts is directly from function stco.findCoasts() or stco.importCoasts
                [north, south, east, west] = coasts

                # should change to contourf ?
                plt.contour(north[2].lon, north[2].lat, north[1].astype(np.bool), alpha=0.5, cmap="Greys")
                plt.contour(south[2].lon, south[2].lat, south[1].astype(np.bool), alpha=0.5, cmap="Greys")
                plt.contour(east[2].lon, east[2].lat, east[1].astype(np.bool), alpha=0.5, cmap="Greys")
                plt.contour(west[2].lon, west[2].lat, west[1].astype(np.bool), alpha=0.5, cmap="Greys")

            else:
                # coasts is a list of Fields from a FieldSet
                for field in coasts:
                    plt.contourf(field.grid.lon, field.grid.lat, field.data[0].astype(np.bool), alpha=0.5, cmap="Greys")

        elif isinstance(coasts, FieldSet):
            # coasts is a FieldSet
            try:
                coasts.F_all
                plt.contourf(coasts.F_all.grid.lon, coasts.F_all.grid.lat, coasts.F_all.data[0].astype(np.bool), alpha=0.5, cmap="Greys")
            except NameError:
                print "plotTrajectories(): cannot find coastfields (F_all) in parameter 'coast'"

        else:
            print "plotTrajectories(): 'coasts' is not a list of coasts or a field"

    if field is not None:
        if field.U.grid == field.V.grid:
            lon_ = field.U.grid.lon
            lat_ = field.U.grid.lat
        else:
            print "plotTrajectories(): grids of field.U and field.V are not the same, trying to continue"
            lon_ = field.U.grid.lon
            lat_ = field.U.grid.lat

        ### index [0] is probably wrong !!!
        color = np.hypot(field.U.data[0], field.V.data[0])
        plt.quiver(lon_, lat_, field.U.data[0], field.V.data[0], color)

    if savefile is not None:
        plt.savefig(savefile)
        print "plotTrajectories(): plot saved as '{}'".format(savefile)
    if show:
        print "plotTrajectories(): showing plot"
        plt.show()
