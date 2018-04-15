""" Functions to make plots for stuck particles

plotfunction: plotVelocity(vector, field="U", coords=None, show=None, savefile=None, vmax=None)
    np.
    math.
    plt.
plotfunction: plotAbsoluteVelocity(vector, coords=None, show=None, savefile=None, vmax=None)
plotfunction: showCoast(fields, type=np.bool, show=None, savefile=None, field="all")
plotfunction: plotLocations(subdata, title="", initial=False, show=None, savefile=None, coastfields=None, coasttype=np.bool)
plotfunction: plotHistogram(subdata, width=1, show=None, savefile=None)
"""
import numpy as np
import math

try:
    import matplotlib.pyplot as plt
except:
    print "ERROR: matplotlib not found"


def plotVelocity(vector, field="U", coords=None, show=None, savefile=None, vmax=None):
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
    plt.plot(plon, plat, 'o', markersize=5, color="red")
    plt.plot(mesh_x.flatten(), mesh_y.flatten(), 'o', markersize=3, color="black")
    plt.title("field {} at t={}".format(field, time))
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


def plotLocations(subdata, title="", initial=False, show=None, savefile=None, coastfields=None, coasttype=np.bool):
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
    n_bins = np.arange(0, np.max(time_stuck+time_moving)+5, 5)

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
