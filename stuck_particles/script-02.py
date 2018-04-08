import numpy as np

from multiprocessing import Pool
from multiprocessing import cpu_count

import GlobalStuckParticles as gsp


def gsp_multi(coast_args):
    coasts = coast_args.tolist()[:-1]
    savename = coast_args.tolist()[-1]
    print "Start:", savename
    simulation = [300, 10] #days, minutes
    particles = [np.linspace(-175, 175, num=11), np.linspace(-75, 75, num=9)]

    data = gsp.GlobalStuckParticlesAdvection(simulation, particles, coasts, savename)

    gsp.GlobalStuckParticlesAnalysis([0, 0], savename, particledata=None, coast=True, plot=True, histogram=True, show=False, locations=True, velocities=False)

    print "Stop:", savename


coast_list = []
coast_list.append(np.array([False, -1, 0, "GlobalStuckParticles_F_-10_00/"], dtype=object))
coast_list.append(np.array([True, -1, 0, "GlobalStuckParticles_T_-10_00/"], dtype=object))
coast_list.append(np.array([True, -2, 0, "GlobalStuckParticles_T_-20_00/"], dtype=object))
coast_list.append(np.array([True, -0.5, 0, "GlobalStuckParticles_T_-05_00/"], dtype=object))
coast_list.append(np.array([True, 0, 1, "GlobalStuckParticles_T_-00_10/"], dtype=object))
coast_list.append(np.array([True, 0, 2, "GlobalStuckParticles_T_-00_20/"], dtype=object))
coast_list.append(np.array([True, 0, 0.5, "GlobalStuckParticles_T_-00_05/"], dtype=object))

pool = Pool(3)
pool.map(gsp_multi, coast_list)
