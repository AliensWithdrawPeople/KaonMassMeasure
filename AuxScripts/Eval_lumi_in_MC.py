import numpy as np

energy = np.array([505, 508, 508.5, 509, 509.5, 510, 510.5, 511, 511.5, 514])
entries = np.array([846439, 745445, 754083, 2350284, 2884456, 2772953, 787669, 783597, 795083, 893345])
cross_section = np.array([32.3232, 300.888, 522.891, 936.339, 1395.78, 1342.44, 952.039, 583.632, 424.024, 119.118])
cross_section = np.array([32.3232, 300.888, 522.891, 936.339, 1395.78, 1342.44, 952.039, 583.632, 424.024, 119.118])

lumi = list(np.round(entries / cross_section, 3))
lumi = [454.9, 549.0, 1395.1, 979.5, 2621.6, 2041.3, 1053.5, 780.3, 561.3, 520.9]
print(lumi)
print(sum(lumi))