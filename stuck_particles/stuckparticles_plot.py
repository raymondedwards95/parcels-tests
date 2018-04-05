""" Functions to make plots for stuck particles

plotfunction: plotVelocity(vector, field="U", coords=None, show=None, savefile=None, vmax=None)
    np.
    math.
    plt.
plotfunction: plotAbsoluteVelocity(vector, coords=None, show=None, savefile=None, vmax=None)
plotfunction: showCoast(fields, type=np.bool, origin="GlobCurrent", show=None, savefile=None)
plotfunction: plotLocations(subdata, title="", initial=False, show=None, savefile=None, coastfield=None, coastorigin="GlobCurrent", coasttype=np.bool)
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
    if show:
        plt.show()
    if savefile is not None:
        plt.savefig(savefile)


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
    if show:
        plt.show()
    if savefile is not None:
        plt.savefig(savefile)


def showCoast(fields, type=np.bool, origin="GlobCurrent", show=None, savefile=None):
    """ Plot values in fields in a contour plot """
    if savefile is None:
        show = True
    if show is None and savefile is not None:
        show = False

    if len(np.shape(fields)) == 3 and np.shape(fields)[0] == 2:
        field = fields[0] + fields[1]
    elif len(np.shape(fields)) == 2:
        field = fields
    else:
        print "showCoast(): fields not in correct form"
        return

    field = field.astype(type)

    ny, nx = np.shape(field)
    if origin == "GlobCurrent":
        lons = np.linspace(-181.125, 181.125, num=nx)
        lats = np.linspace(-79.875, 79.875, num=ny)
    else:
        lons = np.arange(nx)
        lats = np.arange(ny)

    plt.figure()
    plt.contourf(lons, lats, field)
    plt.title("Locations of coasts")
    plt.xlabel("longitude")
    plt.ylabel("latitude")
    plt.colorbar()
    plt.grid()
    if show:
        print "showCoast(): showing plot"
        plt.show()
    if savefile is not None:
        plt.savefig(savefile)


def plotLocations(subdata, title="", initial=False, show=None, savefile=None, coastfield=None, coastorigin="GlobCurrent", coasttype=np.bool):
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

    if coastfield is not None:
        if len(np.shape(coastfield)) == 3 and np.shape(coastfield)[0] == 2:
            field = coastfield[0] + coastfield[1]
        elif len(np.shape(coastfield)) == 2:
            field = coastfield
        else:
            print "plotLocations(): coastfield not in correct form"
            coastfield = None

        field = field.astype(coasttype)

        ny, nx = np.shape(field)
        if coastorigin == "GlobCurrent":
            lons = np.linspace(-181.125, 181.125, num=nx)
            lats = np.linspace(-79.875, 79.875, num=ny)
        else:
            lons = np.arange(nx)
            lats = np.arange(ny)


    colors = ["red", "blue", "green", "pink", "orange", "purple"]
    number = len(colors)

    plt.figure()
    for i in range(len(subdata)):
        m = i%number

        plt.plot(subdata[i][1], subdata[i][2], "o", markersize=3, color=colors[m], label="Particle {} at ({:.1f}, {:.1f})".format(int(subdata[i][0]), subdata[i][1], subdata[i][2]))

        if initial and subdata[0][-1] > 2:
            plt.plot(subdata[i][5], subdata[i][6], "o", markersize=1, color=colors[m])
            plt.plot([subdata[i][1], subdata[i][5]], [subdata[i][2], subdata[i][6]], "--", color=colors[m])

    if coastfield is not None:
        plt.contourf(lons, lats, field, alpha=0.5, cmap="Greys")
        # plt.colorbar()

    plt.grid()
    plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    plt.xlabel("longitude")
    plt.ylabel("latitude")
    plt.title(title)

    if show:
        plt.show()
    if savefile is not None:
        plt.savefig(savefile)


def plotHistogram(subdata, width=1, show=None, savefile=None):
    """ Show number of stuck particles in a histogram using data created in
    stuckparticles_analysis.
    """
    if savefile is None:
        show = True
    if show is None and savefile is not None:
        show = False

    if subdata[0][-1] < 2:
        print "plotHistogram(): missing grid information."
        return

    list = []
    for p in subdata:
        list.append(p[3]/(24.*60.*60.))

    n_bins = int(round((np.max(list) - np.min(list)) / width))

    plt.figure()
    plt.hist(list, n_bins, facecolor="green", alpha=0.75)
    plt.xlabel("days stuck")
    plt.ylabel("number of particles")
    plt.title("")
    plt.grid()

    if show:
        print "plotHistogram(): showing plot"
        plt.show()
    if savefile is not None:
        plt.savefig(savefile)