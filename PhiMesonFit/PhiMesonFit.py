import numpy as np
from xsection import Xsection
import uproot as up
import scipy

E = np.array([501, 503, 505, 508, 508.5, 509, 509.5, 510, 510.5, 511, 511.5, 514, 517, 520, 525, 530])
energy_points = ["501", "503", "505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514", "517", "520", "525", "530"]
energy = np.array([500.9293, 502.9799, 504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 511.4479, 513.864, 516.9, 519.897, 524.902, 529.999])
# lumi = np.array([519.365, 1232.43, 454.856, 548.972, 1395.09, 979.506, 2621.62, 2041.31, 1053.46, 780.307, 561.27, 520.924])
# lumi_err = np.array([0.996744, 1.54306, 0.94111, 1.0394, 1.6559, 1.38908, 2.27821, 2.01353, 1.44786, 1.24862, 1.06109, 1.02696])

lumi = np.array([463.81, 1213.94, 310.388, 512.997, 1364.79, 962.61, 2512.28, 1963.85, 1039.19, 697.626, 509.562, 505.252, 533.464, 481.78, 583.954, 537.448])
lumi_err = np.array([1.05125, 1.70613, 0.866127, 1.12028, 1.8298, 1.53774, 2.48812, 2.20256, 1.60331, 1.31385, 1.125, 1.12422, 1.16215, 1.11148, 1.23467, 1.19585])

                
e_MC = np.array([0.26591, 0.2577, 0.2517, 0.24325, 0.24028, 0.23805, 0.23718, 0.23539, 0.23547, 0.23485, 0.23438, 0.23102, 0.22891, 0.22683, 0.22242, 0.21371])
# e_MC = np.array([0.26591, 0.2577, 0.2517, 0.24325, 0.24028, 0.23805, 0.23718, 0.23539, 0.23547, 0.23485, 0.23438, 0.22889])
# e_MC = np.array([0.25902, 0.25123, 0.24398, 0.24325, 0.24028, 0.23805, 0.23718, 0.23539, 0.23547, 0.23485, 0.23438, 0.22889])
e_MC_err = np.sqrt(e_MC * (1 - e_MC) * 1e-6)
n_events = [302, 2620, 1735, 27365, 124438, 151469, 591058, 456543, 191690, 85594, 48069, 16469, 8922, 5225, 4018, 2559]

# n_events = []
# for point in energy_points:
#     with up.open(f"C:/work/Science/BINP/Kaon Mass Measure/tr_ph/PhiXSection/Kn{point}.root:Kn") as KnTree: # type: ignore
#         n_events.append(len(KnTree['emeas'].array())) # type: ignore
        
efficiency = e_MC

xsec_vis = n_events / lumi

xsec_err_energy_part = np.array([0.00903471, 0.015707, 0.0667149, 1.00209, 6.51135, 4.1341, 3.67391, 6.89706, 4.01397, 1.47482, 1.1773, 0.204747, 0.137859, 0.0518114, 0.019589, 0.0228199])
xsec_vis_err = np.sqrt( (np.sqrt(n_events) / lumi / efficiency)**2 + (n_events / efficiency / (lumi**2))**2 * lumi_err**2  * (n_events / lumi / (efficiency**2))**2 * e_MC_err**2  + xsec_err_energy_part**2)
xsec = n_events / lumi / efficiency

print("E_beam =", energy_points)
print("N_events =", n_events)
print("eff =", list(np.round(efficiency, 6)))
print("eff_err =", list(np.round(e_MC_err, 6)))
print("sigma_vis =", list(np.round(n_events / lumi, 5)))
print("sigma_vis (efficiency accounted)=", list(np.round(xsec, 5)))
print("sigma^(err)_vis =", list(np.round( xsec_vis_err, 5)) )

# vis_xsection = Xsection(np.array(E), xsec, (500, 515))
# a, b = vis_xsection.Get_Data()