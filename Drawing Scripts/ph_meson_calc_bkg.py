import numpy as np

e = np.array([1010.47, 1012.96, 1015.07, 1016.11, 1017.15, 1017.16, 1018.05, 1019.12, 1019.21, 1019.40, 1019.90, 1021.22, 1021.31, 1022.08, 1022.74, 1023.26, 1025.32, 1027.96, 1029.09, 1033.91 , 1040.03, 1049.86, 1050.86, 1059.95])
k_ch_xsec = np.array([69.87 * 0.735, 152.45 * 0.728, 341.10 * 0.718, 575.08 * 0.712, 993.19 * 0.706, 984.71 * 0.706, 1584.27 * 0.706, 2228.59 * 0.721, 2230.81 * 0.724, 2233.66 * 0.730, 2127.07 * 0.752, 1325.01 * 0.829, 1308.31 * 0.835, 933.95 * 0.885, 710.23 * 0.928, 595.03 * 0.961, 334.77 * 1.077, 191.64 * 1.200, 159.94 * 1.244,89.65 * 1.392, 55.87 * 1.509, 34.47 * 1.604, 33.89 * 1.609, 25.93 * 1.640])
print(k_ch_xsec)
k_ch_bkg = k_ch_xsec * 0.000256
pi_neutral1 = np.array([6.32] * len(e[e<1020])) * 0.04
pi_neutral2 = np.array([8.09] * len(e[(e > 1020) & (e < 1040)])) * 0.04
pi_neutral3 = np.array([9.85] * len(e[e > 1040])) * 0.04
pi_neutral = np.concatenate((pi_neutral1, pi_neutral2, pi_neutral3))

pi_charged1 = np.array([1.2] * len(e[e < 1020])) * 0.07
pi_charged2 = np.array([1.5] * len(e[(e > 1020) & (e < 1040)])) * 0.07
pi_charged3 = np.array([2] * len(e[e > 1040])) * 0.07
pi_charged = np.concatenate((pi_charged1, pi_charged2, pi_charged3))

print(list(np.round(e * 1e-3, 4)))
print(list(np.round(k_ch_bkg + pi_neutral + pi_charged, 3)))
print(list(np.round(pi_neutral, 3)))

