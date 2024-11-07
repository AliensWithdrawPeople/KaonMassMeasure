import numpy as np
from xsection import Xsection
import uproot as up
import pandas as pd

df = pd.read_csv('PhiMesonFit/track_corr.dat', sep=' ')
cs_track_rec_corr_dict = dict(zip([el.split('_')[0] for el in df['elabel']], 
                             zip(df['full_cs_corr'], df['full_cs_corr_er'])
                             )
                         )
cs_track_rec_corr = np.array(df['full_cs_corr'])
cs_track_rec_corr_err = np.array(df['full_cs_corr_er'])
print(list(cs_track_rec_corr))
print(list(cs_track_rec_corr_err))

E = np.array([501, 503, 505, 508, 508.5, 509, 509.5, 510, 510.5, 511, 511.5, 514, 517, 520, 525, 530])
energy_points = ["501", "503", "505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514", "517", "520", "525", "530"]
energy = np.array([500.9293, 502.9799, 504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 511.4479, 513.864, 516.9, 519.897, 524.902, 529.999])

lumi_v9 = np.array([518.56, 1225.95, 338.713, 548.328, 1378.21, 975.715, 2611.47, 2022.2, 1051.8, 704.03, 560.555, 509.673, 537.708, 488.13, 600.65, 543.963])
lumi_v9_err = np.array([0.996242, 1.53986, 0.812959, 1.03917, 1.64823, 1.38669	, 2.27446, 2.00449, 1.44699, 1.18695, 1.06076, 1.01584, 1.05068, 1.00447, 1.1248, 1.08156])

lumi_v9_by_runs = np.array([468.021, 1225.95, 311.608, 515.718, 1374.09, 970.723, 2536.03, 1986.74, 1051.8, 704.028, 514.846, 509.672, 537.708, 488.13, 591.476, 543.961])
lumi_v9_by_runs_err = np.array([0.946451, 1.53986, 0.779753, 1.00779, 1.64577, 1.38314, 2.24137, 1.98684, 1.44699, 1.18695, 1.01659, 1.01584, 1.05068, 1.00447, 1.11618, 1.08156])

lumi_v9_by_runs_new = np.array([463.81, 1213.94, 310.388, 512.997, 1364.79, 962.61, 2512.28, 1963.85, 1039.19, 697.626, 509.562, 505.252, 533.464, 481.78, 583.954, 537.448])
lumi_v9_by_runs_new_err = np.array([1.05125, 1.70613, 0.866127, 1.12028, 1.8298, 1.53774, 2.48812, 2.20256, 1.60331, 1.31385, 1.125, 1.12422, 1.16215, 1.11148, 1.23467, 1.19585])

lumi = lumi_v9_by_runs_new           
lumi_err = lumi_v9_by_runs_new_err
           
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
# New: Const bkg (only right sideband); For E >= 517 sideband 4pi
n_events = [260, 2579, 1781, 27053, 123149, 150554, 587751, 453786, 190345, 84830, 47554, 17018.021, 9184.62, 5388.09, 3967.18, 2498.96]
            
n_events = [267, 2569, 1781, 26862, 121999, 148571, 580460, 448154, 187365, 84091, 47117, 17018.021, 9184.62, 5388.09, 3967.18, 2498.96]

events_total = []
# for i, point in enumerate(energy_points):
#     with up.open(f"C:/work/Science/BINP/Kaon Mass Measure/tr_ph/PhiXSection/Kn{point}.root:Kn") as KnTree: # type: ignore
#         events_total.append(len(KnTree['emeas'].array())) # type: ignore
# events_total = np.array(events_total)

# New: MC + pol2 fit 
# ev_sig = np.array([303.308, 2664.73, 1792.23, 27217.7, 124000, 151043, 589622, 455541, 191265, 85200.7, 47883.5, 17112.1, 9295.12, 5465.21, 4224.31, 2660.99])
# ev_sig_err = np.array([18.7858, 53.0068, 43.1968, 166.921, 352.878, 52.2047, 3.90865e-05, 68.892, 440.145, 278.756, 7.9838e-07, 106.129, 98.6353, 76.1091, 67.4571, 53.9637])
# ev_bkg = np.array([281.694, 492.196, 122.718, 120.75, 329.694, 210.224, 745.039, 403.48, 135.046, 249.214, 215.988, 287.243, 398.87, 376.729, 573.031, 573.982])
# ev_bkg_err = np.array([18.2023, 25.2735, 14.0677, 27.7477, 29.0511, 94.7366, 3.39091e-06, 7.3185, 50.7135, 31.2297, 1.30522e-07, 24.2181, 28.8786, 26.6231, 30.2251, 28.7813])
ev_sig = np.array([307.566, 2654.43, 1797.53, 27817.9, 126417, 154515, 600473, 465934, 194871, 86986.6, 47916.4, 17060.1, 9309.34, 5724.98, 4172.94, 2689.46])
ev_sig_err = np.array([18.9421, 53.0088, 41.8419, 166.802, 355.696, 397.511, 777.746, 690.166, 446.172, 294.985, 213.353, 127.298, 95.248, 79.2939, 64.8708, 55.4257])
sig_extraction_syst = np.array([0.022, 0.009, 0.002, 0.03, 0.001, 0.007, 0.019, 0.008, 0.006, 0.001, 0.011, 0.001, 0.002, 0.011, 0.014, 0.016])
ev_sig_err = ev_sig * ((ev_sig_err / ev_sig)**2 + sig_extraction_syst**2)**0.5

ev_bkg = np.array([297.766, 537.771, 111.102, 239.959, 1106.63, 724.437, 2439.95, 2035.55, 973.33, 617.743, 386.658, 276.97887, 397.379, 309.905, 410.821, 320.037])
ev_bkg_err = ev_bkg**0.5

xsec_bkg: np.ndarray = np.round((ev_bkg / lumi), 5)
xsec_bkg_err = np.round(np.sqrt( (1 / lumi)**2 * ev_bkg_err**2 + (ev_bkg / (lumi**2))**2 * lumi_err**2), 5)

bkg_events_estimated = lumi * np.array([0.346154, 0.37076, 0.3926, 0.474208, 0.538273, 0.659605, 0.791, 0.792818, 0.693963, 0.59955, 0.546487, 0.447324, 0.417396, 0.407344, 0.399081, 0.396])
xsec_bkg_estimated = np.array([0.346154, 0.37076, 0.3926, 0.474208, 0.538273, 0.659605, 0.791, 0.792818, 0.693963, 0.59955, 0.546487, 0.447324, 0.417396, 0.407344, 0.399081, 0.396])
# print(list(np.round(xsec_bkg, 3)))
# print(list(np.round(xsec_bkg_estimated, 3)))
# print(list(np.round(100 * (xsec_bkg_estimated) / xsec_bkg, 2)))
# n_events = events_total - bkg_events_estimated
n_events = ev_sig
n_events_err = ev_sig_err
efficiency = e_MC

xsec_vis = n_events / lumi

xsec_err_energy_part = np.array([0.00903471, 0.015707, 0.0667149, 1.00209, 6.51135, 4.1341, 3.67391, 6.89706, 4.01397, 1.47482, 1.1773, 0.204747, 0.137859, 0.0518114, 0.019589, 0.0228199])
xsec_vis_err = np.sqrt( (n_events_err / lumi / efficiency)**2 + 
                       (n_events / efficiency / (lumi**2))**2 * lumi_err**2  + 
                       (n_events / lumi / (efficiency**2))**2 * e_MC_err**2 +
                       0 * (n_events / efficiency / lumi / (cs_track_rec_corr**2))**2 * cs_track_rec_corr_err**2
                        )
# xsec = n_events / lumi / efficiency / cs_track_rec_corr
xsec = n_events / lumi / efficiency

print("E_beam =", energy_points)
# print("N_events =", list(n_events))
# print("eff =", list(np.round(efficiency, 6)))
# print("eff_err =", list(np.round(e_MC_err, 6)))
# print("n_MC =", list(n_MC))
# print("sigma_vis =", list(np.round(n_events / lumi, 5)))
print("sigma_vis (efficiency accounted)=", list(np.round(xsec, 5)))
print("sigma^(err)_vis =", list(np.round( xsec_vis_err, 6)) )
# print("xsec_bkg =", list(xsec_bkg))
# print("xsec_bkg_err =", list(xsec_bkg_err))

# vis_xsection = Xsection(np.array(E), xsec, (500, 515))
# a, b = vis_xsection.Get_Data()




energy = np.array([500.9293, 502.9799, 504.8, 507.862, 508.404, 508.957, 509.528, 509.956, 510.458, 511.035, 511.4479, 513.864, 516.9, 519.897, 524.902, 529.999])

ev_sig_m1 = np.array([260, 2579, 1781, 27053, 123149, 150554, 587751, 453786, 190345, 84830, 47554, 16892, 9138, 5349, 4059, 2555])
ev_sig_m1_err = ev_sig_m1**0.5
ev_bkg_m1 = np.array([297.766, 537.771, 111.102, 239.959, 1106.63, 724.437, 2439.95, 2035.55, 973.33, 617.743, 386.658, 276.97887, 397.379, 309.905, 410.821, 320.037])
ev_bkg_m1_err = ev_bkg**0.5

ev_sig_m2 = np.array([307.566, 2654.43, 1797.53, 27817.9, 126417, 154515, 600473, 465934, 194871, 86986.6, 47916.4, 17060.1, 9309.34, 5724.98, 4172.94, 2689.46])
ev_sig_m2_err = np.array([18.9421, 53.0088, 41.8419, 166.802, 355.696, 397.511, 777.746, 690.166, 446.172, 294.985, 213.353, 127.298, 95.248, 79.2939, 64.8708, 55.4257])
ev_bkg_m2 = np.array([281.694, 492.196, 122.718, 120.75, 329.694, 210.224, 745.039, 403.48, 135.046, 249.214, 215.988, 287.243, 398.87, 376.729, 573.031, 573.982])
ev_bkg_m2_err = np.array([18.2023, 25.2735, 14.0677, 27.7477, 29.0511, 94.7366, 3.39091e-06, 7.3185, 50.7135, 31.2297, 1.30522e-07, 24.2181, 28.8786, 26.6231, 30.2251, 28.7813])


df = pd.DataFrame({ 'Energy, MeV' : energy, 
                    '$N_{events}$ for #1': ev_sig_m1,
                    r'$\delta N_{events}$ for #1': ev_sig_m1_err,
                    '$N_{bkg}$ for #1' : ev_bkg_m1,
                    r'$\delta N_{bkg}$ for #1': ev_bkg_m1_err,
                    '$N_{events}$ for #2': ev_sig_m2,
                    r'$\delta N_{events}$ for #2': ev_sig_m2_err,
                    '$N_{bkg}$ for #2' : ev_bkg_m2,
                    r'$\delta N_{bkg}$ for #2': ev_bkg_m2_err
                    })

df.to_latex('tmp.txt', float_format="%.3f", longtable=False, index=False)