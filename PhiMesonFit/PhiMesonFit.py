import numpy as np
from xsection import Xsection
import uproot as up
import scipy

E = np.array([501, 503, 505, 508, 508.5, 509, 509.5, 510, 510.5, 511, 511.5, 514, 517, 520, 525, 530])
energy_points = ["501", "503", "505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514", "517", "520", "525", "530"]
energy = np.array([500.9293, 502.9799, 504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 511.4479, 513.864, 516.9, 519.897, 524.902, 529.999])
# lumi_v8 = np.array([519.365, 1232.43, 454.856, 548.972, 1395.09, 979.506, 2621.62, 2041.31, 1053.46, 780.307, 561.27, 520.924])
# lumi_v8_err = np.array([0.996744, 1.54306, 0.94111, 1.0394, 1.6559, 1.38908, 2.27821, 2.01353, 1.44786, 1.24862, 1.06109, 1.02696])

lumi_v8_by_runs = np.array([463.81, 1213.94, 310.388, 512.997, 1364.79, 962.61, 2512.28, 1963.85, 1039.19, 697.626, 509.562, 505.252, 533.464, 481.78, 583.954, 537.448])
lumi_v8_by_runs_err = np.array([1.05125, 1.70613, 0.866127, 1.12028, 1.8298, 1.53774, 2.48812, 2.20256, 1.60331, 1.31385, 1.125, 1.12422, 1.16215, 1.11148, 1.23467, 1.19585])

lumi_v8 = np.array([513.732, 1213.94, 337.389, 545.408, 1368.86, 967.534, 2586.77, 1999.12, 1039.19, 697.626, 554.614, 505.252, 533.464, 481.78, 593.037, 537.448])
lumi_v8_err = np.array([1.10628, 1.70613, 0.903094, 1.15518, 1.83251, 1.54165, 2.52485, 2.22226, 1.60331, 1.31385, 1.17358, 1.12422, 1.16215, 1.11148, 1.24425, 1.19585])

lumi_v9 = np.array([518.56, 1225.95, 338.713, 548.328, 1378.21, 975.715, 2611.47, 2022.2, 1051.8, 704.03, 560.555, 509.673, 537.708, 488.13, 600.65, 543.963])
lumi_v9_err = np.array([0.996242, 1.53986, 0.812959, 1.03917, 1.64823, 1.38669	, 2.27446, 2.00449, 1.44699, 1.18695, 1.06076, 1.01584, 1.05068, 1.00447, 1.1248, 1.08156])

lumi_v9_by_runs = np.array([468.021, 1225.95, 311.608, 515.718, 1374.09, 970.723, 2536.03, 1986.74, 1051.8, 704.028, 514.846, 509.672, 537.708, 488.13, 591.476, 543.961])
lumi_v9_by_runs_err = np.array([0.946451, 1.53986, 0.779753, 1.00779, 1.64577, 1.38314, 2.24137, 1.98684, 1.44699, 1.18695, 1.01659, 1.01584, 1.05068, 1.00447, 1.11618, 1.08156])

lumi = lumi_v9_by_runs           
lumi_err = lumi_v9_by_runs_err
           
e_MC = np.array([0.25922, 0.25138, 0.24396, 0.23382, 0.23255, 0.23099, 0.23027, 0.22983, 0.22779, 0.22737, 0.22637, 0.22331, 0.21986, 0.21707, 0.21214, 0.20391])
e_MC = np.array([0.25945, 0.25114, 0.2444, 0.2339, 0.23235, 0.23133, 0.23111, 0.22996, 0.22806, 0.22772, 0.22709, 0.22395, 0.22032, 0.21754, 0.21277, 0.2047])
n_MC = []
for point in energy_points:
    with up.open(f"C:/work/Science/BINP/Kaon Mass Measure/tr_ph/mcgpj/tr_ph v9 new form factor/Merged/xsec_cutted/MCGPJ_kskl{point}_Merged_XsecConv.root:tr_ph_merged") as tr_ph_merged: # type: ignore
        n_MC.append(len(tr_ph_merged['emeas'].array())) # type: ignore
e_MC_err = np.sqrt(e_MC * (1 - e_MC) / np.array(n_MC))


# Old: Const bkg
n_events = [305, 2624, 1736, 26421, 120322, 146898, 573681, 442909, 185693, 82720, 46215, 16475, 8931, 5228, 4045, 2550]
# New: Const bkg
n_events = [309, 2630, 1738, 26458, 120514, 147107, 574569, 443659, 186014, 82867, 46320, 16524, 8963, 5252, 4062, 2571]
# New: Const bkg (only right sideband)
n_events = [260, 2579, 1781, 27053, 123149, 150554, 587751, 453786, 190345, 84830, 47554, 16892, 9138, 5349, 4059, 2555]
# for point in energy_points:
#     with up.open(f"C:/work/Science/BINP/Kaon Mass Measure/tr_ph/PhiXSection/Kn{point}.root:Kn") as KnTree: # type: ignore
#         n_events.append(len(KnTree['emeas'].array())) # type: ignore

# New: MC + pol2 fit 
ev_sig = np.array([303.308, 2664.73, 1792.23, 27217.7, 124000, 151043, 589622, 455541, 191265, 85200.7, 47883.5, 17112.1, 9295.12, 5465.21, 4224.31, 2660.99])
ev_sig_err = np.array([18.7858, 53.0068, 43.1968, 166.921, 352.878, 52.2047, 3.90865e-05, 68.892, 440.145, 278.756, 7.9838e-07, 106.129, 98.6353, 76.1091, 67.4571, 53.9637])
ev_bkg = np.array([281.694, 492.196, 122.718, 120.75, 329.694, 210.224, 745.039, 403.48, 135.046, 249.214, 215.988, 287.243, 398.87, 376.729, 573.031, 573.982])
ev_bkg_err = np.array([18.2023, 25.2735, 14.0677, 27.7477, 29.0511, 94.7366, 3.39091e-06, 7.3185, 50.7135, 31.2297, 1.30522e-07, 24.2181, 28.8786, 26.6231, 30.2251, 28.7813])
xsec_bkg: np.ndarray = np.round((ev_bkg / lumi), 5)
xsec_bkg_err = np.round(np.sqrt( (1 / lumi)**2 * ev_bkg_err**2 + (ev_bkg / (lumi**2))**2 * lumi_err**2), 5)

n_events = ev_sig
efficiency = e_MC

xsec_vis = n_events / lumi

xsec_err_energy_part = np.array([0.00903471, 0.015707, 0.0667149, 1.00209, 6.51135, 4.1341, 3.67391, 6.89706, 4.01397, 1.47482, 1.1773, 0.204747, 0.137859, 0.0518114, 0.019589, 0.0228199])
xsec_vis_err = np.sqrt( (np.sqrt(n_events) / lumi / efficiency)**2 + (n_events / efficiency / (lumi**2))**2 * lumi_err**2  + (n_events / lumi / (efficiency**2))**2 * e_MC_err**2)
xsec = n_events / lumi / efficiency / (1-0.0025) / (1-0.0015)

print("E_beam =", energy_points)
print("N_events =", list(n_events))
print("eff =", list(np.round(efficiency, 6)))
print("eff_err =", list(np.round(e_MC_err, 6)))
print("n_MC =", list(n_MC))
print("sigma_vis =", list(np.round(n_events / lumi, 5)))
print("sigma_vis (efficiency accounted)=", list(np.round(xsec, 5)))
print("sigma^(err)_vis =", list(np.round( xsec_vis_err, 6)) )
print("xsec_bkg =", list(xsec_bkg))
print("xsec_bkg_err =", list(xsec_bkg_err))

# vis_xsection = Xsection(np.array(E), xsec, (500, 515))
# a, b = vis_xsection.Get_Data()