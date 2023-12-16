import subprocess as sub
import uproot as up
import numpy as np
import re
from functools import reduce

energy_points = ["501", "503", "505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514", "517", "520", "525", "530"]
# energy_points = ["508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514"]
energy_points = ["509.5"]
aux1, aux2 = "\"", "\\"
eff = []
eff_err = []
pattern = r'n_events = (\d+)'
pattern_bckg = r'N_bckg = ([\d.]+)\r\nN_bckg_err = ([\d.]+)'

n_events = dict()
n_bckg = []
n_bckg_err = []

for en in energy_points:
    command = f"root -l -q \"C:/work/Science/BINP/Kaon Mass Measure/up to date scripts for cutting/cutters/cutter_MC.cpp(\\{aux1 + en + aux2}\")\""
    res = sub.run(command, capture_output=True)
    output = res.stderr if res.stderr else res.stdout    
    print(res.stdout.decode()[118:])
    print("ERROR!", res.stderr.decode())
    
    match = re.search(pattern, res.stdout.decode())
    if match:
        events = int(match.group(1))
        n_events[en] = events

for en in energy_points:
    print(en)
    N_event1 = n_events[en]
    with up.open(f"C:/work/Science/BINP/Kaon Mass Measure/tr_ph/mcgpj/tr_ph v9 new form factor/Merged/xsec_cutted/MCGPJ_kskl{en}_Merged_XsecConv.root:tr_ph_merged") as ksTree: # type: ignore
        tree = ksTree.arrays(['emeas'], library="np") # type: ignore
    N_event2 = len(tree['emeas'])
    eff.append(N_event1 / N_event2)    
    eff_err.append((eff[-1] * (1 - eff[-1])/N_event2)**0.5)    
print("eff =", list(np.round(eff, 5)))
print("eff_err =", list(np.round(eff_err, 5)))

print("n_bckg =", n_bckg)
print("n_bckg_err =", n_bckg_err)
print("N_background total =", reduce(lambda x, y: x + y, n_bckg, 0))
