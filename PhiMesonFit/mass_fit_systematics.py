import numpy as np
import pandas as pd

E = [501, 503, 505, 508, 508.5, 509, 509.5, 510, 510.5, 511, 511.5, 514, 517, 520, 525, 530]

# Base:
ev_sig = np.array([307.566, 2654.43, 1797.53, 27817.9, 126417, 154515, 600473, 465934, 194871, 86986.6, 47916.4, 17060.1, 9309.34, 5724.98, 4172.94, 2689.46])
ev_sig_err = np.array([18.9421, 53.0088, 41.8419, 166.802, 355.696, 397.511, 777.746, 690.166, 446.172, 294.985, 213.353, 127.298, 95.248, 79.2939, 64.8708, 55.4257])

# Bkg variation (pol0):
pol0_ev_sig = np.array([307.865, 2650.92, 1783.34, 26995.5, 121782, 149089, 582332, 449784, 188197, 84428.8, 47331.1, 17001, 9226.62, 5388.75, 4165.12, 2643.74])
pol0_ev_sig_err = np.array([16.505, 52.8282, 25.6167, 165.44, 32.6951, 388.081, 766.773, 673.66, 435.947, 3.21161, 219.393, 132.432, 98.2266, 75.6543, 67.1748, 54.0948])

# Bkg variation (pol1):
pol1_ev_sig = np.array([310.161, 2677.44, 1794.3, 27812.6, 126389, 154164, 589052, 465747, 194850, 86966.3, 47926.9, 17057.4, 9298.27, 5786.4, 4167.36, 2677.65])
pol1_ev_sig_err = np.array([19.0858, 51.9339, 43.2142, 166.759, 355.584, 395.717, 746.566, 689.629, 446.105, 294.891, 213.632, 127.382, 95.1587, 79.3883, 64.8326, 55.1667])

# Bkg variation (pol3):
pol3_ev_sig = np.array([314.583, 2653.7, 1794.26, 27818.1, 126384, 153471, 598885, 462447, 193728, 86979.3, 48451.5, 17063.8, 9329.73, 5732.14, 4234.12, 2670.59])
pol3_ev_sig_err = np.array([19.2857, 53.0436, 41.8555, 166.784, 355.486, 390.632, 771.76, 676.653, 439.284, 294.914, 228.94, 127.203, 0.541575, 0.233404, 65.9276, 54.6682])

# Bkg variation (pol4):
pol4_ev_sig = np.array([299.088, 2661.48, 1781.19, 27005, 122834, 149153, 581701, 449885, 187408, 84484.8, 47292.4, 17056.4, 9225.89, 5384.93, 4322.59, 2660.58])
pol4_ev_sig_err = np.array([8.09392, 37.5164, 43.0881, 28.4587, 352.326, 123.333, 769.298, 642.248, 434.543, 114.852, 220.471, 72.9716, 98.411, 75.8917, 0.0494144, 58.7043])

# Range variation (short = (430 - 560)):
short_ev_sig = np.array([309.668, 2655.5, 1795.77, 27016.1, 126512, 154276, 590919, 465027, 194557, 87089.3, 47646, 17076.7, 9291.36, 5747.55, 4182.02, 2647.74])
short_ev_sig_err = np.array([18.9276, 51.5755, 43.3349, 164.852, 356.238, 397.122, 751.738, 687.266, 445.501, 295.722, 210.616, 127.847, 95.1261, 86.578, 67.0577, 54.6337])

# Range variation (super short = (460 - 530)):
super_short_ev_sig = np.array([310.259, 2664.55, 2304.12, 27066.8, 156428, 149270, 582961, 450394, 188430, 84538.7, 47412.7, 17060.6, 9234.99, 5367.71, 5169.37, 3243.39, ])
super_short_ev_sig_err = np.array([19.1046, 53.453, 56.3843, 165.356, 447.943, 250.429, 76.677, 673.608, 398.113, 9.56341, 218.656, 135.241, 99.9356, 76.997, 85.4356, 66.72, ])

df = pd.DataFrame({'E' : E})
def eval_diff(base: np.ndarray, comp: np.ndarray)->list:
    return list(np.round(np.abs(base - comp)/comp * 100, 2))

def eval_diff_err(base: np.ndarray, comp: np.ndarray)->list:
    return list(np.round((base**2 + comp**2)**0.5, 3))

df["Bkg var (pol0), %"] = eval_diff(ev_sig, pol0_ev_sig)
df["Bkg var (pol0) err, %"] = eval_diff_err(ev_sig_err/ev_sig, pol0_ev_sig_err/pol0_ev_sig)

df["Bkg var (pol1), %"] = eval_diff(ev_sig, pol1_ev_sig)
df["Bkg var (pol1) err, %"] = eval_diff_err(ev_sig_err/ev_sig, pol1_ev_sig_err/pol1_ev_sig)

df["Bkg var (pol3), %"] = eval_diff(ev_sig, pol3_ev_sig)
df["Bkg var (pol3) err, %"] = eval_diff_err(ev_sig_err/ev_sig, pol3_ev_sig_err/pol3_ev_sig)

df["Bkg var (pol4), %"] = eval_diff(ev_sig, pol4_ev_sig)
df["Bkg var (pol4) err, %"] = eval_diff_err(ev_sig_err/ev_sig, pol4_ev_sig_err/pol4_ev_sig)

df["Range var (short), %"] = eval_diff(ev_sig, short_ev_sig)
df["Range var (short) err, %"] = eval_diff_err(ev_sig_err/ev_sig, short_ev_sig_err/short_ev_sig)

df["Range var (super short), %"] = eval_diff(ev_sig, super_short_ev_sig)
df["Range var (super short) err, %"] = eval_diff_err(ev_sig_err/ev_sig, super_short_ev_sig_err/super_short_ev_sig)

df['systematics, %'] = np.maximum.reduce([df["Bkg var (pol1), %"], 
                                          df["Bkg var (pol3), %"], 
                                          df["Range var (short), %"]])
print('E_{beam}, MeV =', df['E'].to_list())
print('systematics, % =', df['systematics, %'].to_list())

print('err =', list(np.round(ev_sig * ((ev_sig**0.5/ev_sig)**2 + (df['systematics, %'] * 1e-2)**2)**0.5, 3)))
print('rel err =', list(np.round(((ev_sig**0.5/ev_sig)**2 + (df['systematics, %'] * 1e-2)**2)**0.5, 3)))
print('rel err =', list(np.round(df['systematics, %'] * 1e-2, 3)))


