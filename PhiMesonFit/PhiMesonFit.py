import numpy as np
import uproot as up
from scipy import integrate

from xsection import Xsection


points = ["505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "514"]
efficiency = np.array([0.130427, 0.13764, 0.138855, 0.140421, 0.143378, 0.142905, 0.142406, 0.145685, 0.146109])
B_Kpipi = 0.692
efficiency = efficiency * B_Kpipi
lumi = np.array([454.856, 548.972, 1395.09, 979.506, 2621.62, 2041.31, 1053.46, 780.307, 520.924])
mean_energies = [504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 513.864]
energy_info = dict(zip(points, mean_energies))

N_events = []
for energy_point in points:
    with up.open("C://work/Science/BINP/Kaon Mass Measure/tr_ph/expKsKl/exp" + energy_point +"_v9.root:ksTree") as ksTree:
        tree = ksTree.arrays(['emeas', 'demeas', 'runnum', 'ksdpsi', 'Y'], library="np")
    N_events.append(len(tree['emeas']))
N_events = np.array(N_events)

xsec = np.round(N_events / efficiency / lumi, 5)

vis_xsection = Xsection(np.array(mean_energies), xsec, (500, 515))
print(vis_xsection.Get_xsection(507.5))

vis_xsection.draw()
