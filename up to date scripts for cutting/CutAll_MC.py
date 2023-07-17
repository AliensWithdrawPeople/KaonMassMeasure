import subprocess as sub
import uproot as up
import numpy as np
import re
from functools import reduce

energy_points = ["501", "503", "505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514", "517", "520", "525", "530"]
# energy_points = ["508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514"]
# energy_points = ["510", "514"]
aux1, aux2 = "\"", "\\"
eff = []
pattern = r'n_events = (\d+)'
pattern_bckg = r'N_bckg = ([\d.]+)\r\nN_bckg_err = ([\d.]+)'

n_events = dict()
n_bckg = []
n_bckg_err = []

for en in energy_points:
    command = f"root -l -q \"C:/work/Science/BINP/Kaon Mass Measure/up to date scripts for cutting/cutter_MC.cpp(\\{aux1 + en + aux2}\")\""
    res = sub.run(command, capture_output=True)
    output = res.stderr if res.stderr else res.stdout    
    print(res.stdout.decode()[118:])
    print("ERROR!", res.stderr.decode())
    
    match = re.search(pattern, res.stdout.decode())
    if match:
        events = int(match.group(1))
        n_events[en] = events
    
    match_bckg = re.search(pattern_bckg, res.stdout.decode())
    if match_bckg:
        n_bckg.append(float(match_bckg.group(1)))
        n_bckg_err.append(float(match_bckg.group(2)))

for en in energy_points:
    print(en)
    N_event1 = n_events[en]
    
    with up.open(f"C:/work/Science/BINP/Kaon Mass Measure/tr_ph/mcgpj/tr_ph v9/EnergySmearing/MCGPJ_kskl{en}_Merged.root:tr_ph_merged") as ksTree: # type: ignore
        tree = ksTree.arrays(['emeas'], library="np") # type: ignore
    N_event2 = len(tree['emeas'])
    eff.append(N_event1 / N_event2)    
print("eff =", list(np.round(eff, 5)))

print("n_bckg =", n_bckg)
print("n_bckg_err =", n_bckg_err)
print("N_background total =", reduce(lambda x, y: x + y, n_bckg, 0))