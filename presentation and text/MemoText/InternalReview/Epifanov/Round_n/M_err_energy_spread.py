import numpy as np

E = np.array([505, 508, 508.5, 509, 509.5, 510, 510.5, 511, 514])
sigma_E = np.array([352, 385, 351, 352, 373, 397, 396, 363, 369, 373])
delta_sigma_E = np.array([13, 13, 20, 18, 8, 11, 13, 11, 23, 17])

Espread = sigma_E / 2**0.5
dEspread = delta_sigma_E / 2**0.5

M = np.array(
    [
        497.618,
        497.628,
        497.615,
        497.612,
        497.614,
        497.613,
        497.622,
        497.618,
        497.625,
        497.607,
    ]
)
M_no_spread = np.array(
    [
        497.619,
        497.637,
        497.615,
        497.616,
        497.617,
        497.618,
        497.631,
        497.61,
        497.631,
        497.6,
    ]
)

dM_spread = np.abs(M - M_no_spread)
sigma_M = dM_spread * dEspread / Espread
print(list(np.round(sigma_M, 5)))
# (dM_plus * dEspread / Espread, dM_minus * dEspread / Espread)

dM = np.array([0.004, 0.006, 0.006, 0.004, 0.004, 0.004, 0.007, 0.008, 0.009, 0.016])
M_err = np.array([0.046, 0.018, 0.012, 0.011, 0.01, 0.01, 0.011, 0.015, 0.019, 0.053])
weight = 1 / M_err
dM_avg = np.sum(dM * weight) / np.sum(weight)
print(dM_avg)


M_std = np.array(
    [
        497.618,
        497.625,
        497.614,
        497.612,
        497.615,
        497.612,
        497.618,
        497.616,
        497.621,
        497.615,
    ]
)

M_alt = np.array(
    [
        497.621,
        497.626,
        497.615,
        497.612,
        497.615,
        497.612,
        497.619,
        497.617,
        497.624,
        497.639,
    ]
)

dM = abs(M_std - M_alt)
print(list(np.round(dM * 1e3, 1)))
# dM = [0.004, 0.003, 0.001, 0.003, 0.001, 0.001, 0.001, 0.004, 0.005, 0.007]
dM = [4, 3, 1, 3, 0, 1, 1, 4, 0, 6]
dM_avg = np.sum(dM * weight) / np.sum(weight)
# dM_avg = np.average(dM)

print(dM_avg)
