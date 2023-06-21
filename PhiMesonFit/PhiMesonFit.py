import numpy as np
import uproot as up
from xsection import Xsection
import matplotlib.pyplot as ppl
import matplotlib
matplotlib.use('TkAgg')


points = ["501", "503", "505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514"]
e_MC = np.array([0.100742, 0.14724, 0.168875, 0.189286, 0.191459, 0.194718, 0.19831, 0.198671, 0.199915, 0.202736, 0.202236, 0.205726])
B_Kpipi = 0.692
efficiency = e_MC * B_Kpipi
lumi = np.array([519.365, 1232.43, 454.856, 548.972, 1395.09, 979.506, 2621.62, 2041.31, 1053.46, 780.307, 561.27, 520.924])
lumi_err = np.array([0.996744, 1.54306, 0.94111, 1.0394, 1.6559, 1.38908, 2.27821, 2.01353, 1.44786, 1.24862, 1.06109, 1.02696])
mean_energies = [500.9293, 502.9799, 504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 511.4479, 513.864]
energy_info = dict(zip(points, mean_energies))

N_events = []
for energy_point in points:
    with up.open("C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKsKl/KsKl_Phi/exp" + energy_point +".root:ksTree") as ksTree:
        tree = ksTree.arrays(['emeas', 'demeas', 'runnum', 'ksdpsi', 'Y'], library="np")
    N_events.append(len(tree['emeas']))
N_events = np.array(N_events)

xsec = np.round(N_events / efficiency / lumi, 5)

vis_xsection = Xsection(np.array(mean_energies), xsec, (500, 515))

a, b = vis_xsection.Get_Data()

ppl.plot(a, b)
ppl.show()