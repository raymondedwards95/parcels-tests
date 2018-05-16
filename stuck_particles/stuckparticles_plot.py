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
import numpy as np
import math

from netCDF4 import Dataset

from parcels import Field, FieldSet

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


def plotLocations(subdata, title="", initial=False, show=None, savefile=None, coastfields=None, legend=False, coasttype=np.bool):
    """ Plot locations of particles in subdata or data from stuckparticles_analysis.
    Assuming that the first three values are (id, lon, lat).
    """
    if savefile is None:
        show = True
    if show is None and savefile is not None:
        show = False

    if initial and subdata[0][-1] < 3:
        print "plotLocations(): not enough data."
        initial = False

    if coastfields is not None:
        if len(coastfields) == 4:
            [field_coast_U, field_coast_V, lons, lats] = coastfields
        else:
            print "plotLocations(): 'coastfields' contains not enough information, skipping coasts"
            coastfields = None

        field = np.sqrt(np.square(field_coast_U) + np.square(field_coast_V))
        field = field.astype(coasttype)


    colors = ["red", "blue", "green", "pink", "orange", "purple"]
    number = len(colors)

    plt.figure()
    for i in range(len(subdata)):
        m = i%number

        plt.plot(subdata[i][1], subdata[i][2], "o", markersize=4, color=colors[m], label="Particle {} at ({:.1f}, {:.1f})".format(int(subdata[i][0]), subdata[i][1], subdata[i][2]))

        if initial and subdata[0][-1] > 2:
            plt.plot(subdata[i][5], subdata[i][6], "o", markersize=1, color=colors[m])
            plt.plot([subdata[i][1], subdata[i][5]], [subdata[i][2], subdata[i][6]], "--", color=colors[m])

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


def plotHistogram(subdata, width=1, show=None, savefile=None, title=""):
    """ Show number of stuck particles in a histogram using data created in
    stuckparticles_analysis.
    """
    if savefile is None:
        show = True
    if show is None and savefile is not None:
        show = False

    if subdata[0][-1] < 2:
        print "plotHistogram(): missing time information."
        return

    time_stuck = np.array([])
    time_moving = np.array([])
    for p in subdata:
        time_stuck = np.append(time_stuck, p[3]/(24*60*60))
        time_moving = np.append(time_moving, p[4]/(24*60*60))

    time_stuck_red = np.trim_zeros(sorted(time_stuck, reverse=True))

    # n_bins = int(round((np.max(time_stuck+time_moving) - np.min(time_stuck+time_moving)) / width))
    n_bins = np.arange(0, np.max(time_stuck+time_moving)+width, width)

    plt.figure()
    plt.hist(time_stuck_red, n_bins, facecolor="green", alpha=0.75, edgecolor='black', linewidth=0.1)
    plt.xticks(n_bins)
    plt.xlabel("days stuck")
    plt.ylabel("number of particles")
    plt.title(title+"\n{} particles in histogram".format(len(time_stuck)))
    plt.grid()

    if savefile is not None:
        plt.savefig(savefile)
        print "plotHistogram(): plot saved as '{}'".format(savefile)
    if show:
        print "plotHistogram(): showing plot"
        plt.show()


def scatterStuckMoving(subdata, show=None, savefile=None, title=""):
    """ Create a scatter plot from time_moving and time_stuck in subdata
    Plot shows if a particle was stuck more than once.
    """
    if savefile is None:
        show = True
    if show is None and savefile is not None:
        show = False

    if subdata[0][-1] < 2:
        print "scatterStuckMoving(): missing time information."
        return

    time_stuck = np.array([])
    time_moving = np.array([])
    for p in subdata:
        time_stuck = np.append(time_stuck, p[3]/(24*60*60))
        time_moving = np.append(time_moving, p[4]/(24*60*60))

    x = np.linspace(0, np.max(time_moving+time_stuck), 2)
    y = np.max(time_moving+time_stuck) - x

    plt.figure()
    plt.scatter(time_stuck, time_moving)
    plt.plot(x, y, "--", color="red")
    plt.xlabel("days stuck")
    plt.ylabel("days moved")
    plt.xlim([-0.1, max(time_moving+time_stuck)*1.05])
    plt.ylim([-0.1, max(time_moving+time_stuck)*1.05])
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
