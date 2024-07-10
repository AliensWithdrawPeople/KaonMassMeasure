import subprocess as sub
from functools import reduce
import re
import numpy as np

energy_points = ["501", "503", "505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514", "517", "520", "525", "530"]
# energy_points = ["508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514"]
# energy_points = ["517",  "520", "525", "530"]
aux1, aux2 = "\"", "\\"
pattern = r'luminosity = ([\d.]+)\r\nluminosity_err = ([\d.]+)\r\nentries = (\d+)\r\nrunCounter = (\d+)'
pattern_bckg = r'N_bckg = ([\d.]+)\r\nN_bckg_err = ([\d.]+)'
        
luminosity = []
luminosity_err = []
n_signal = []
run_counter = []

n_bckg = []
n_bckg_err = []

for en in energy_points:
    command = f"root -l -q \"C:/work/Science/BINP/Kaon Mass Measure/up to date scripts for cutting/cutters/cutter.cpp(\\{aux1 + en + aux2}\")\""
    res = sub.run(command, capture_output=True)
    output = res.stderr if res.stderr else res.stdout
    print(res.stdout.decode()[118:])
    
    match = re.search(pattern, res.stdout.decode())
    if match:
        lumi = float(match.group(1))
        lumi_err = float(match.group(2))
        events = int(match.group(3))
        n_signal.append(events)
        luminosity.append(lumi)
        luminosity_err.append(lumi_err)
        run_counter.append(int(match.group(4)))
    
    match_bckg = re.search(pattern_bckg, res.stdout.decode())
    if match_bckg:
        n_bckg.append(float(match_bckg.group(1)))
        n_bckg_err.append(float(match_bckg.group(2)))

print("How many matches?", len(n_signal)) 
print("Is it ok?", "Yes!" if (len(n_signal) == len(energy_points)) else "No...")
print("n_signal =", n_signal)
print("n_bckg =", n_bckg)
print("luminosity =", luminosity)
print("luminosity_err =", luminosity_err)
print("sigma_bckg =", list(np.round(np.array(n_bckg)/np.array(luminosity), 5)))
print("sigma_bckg_err =", list(np.round(np.sqrt((np.sqrt(np.array(n_bckg)) / np.array(luminosity))**2 + (np.array(n_bckg)/np.array(luminosity)**2)**2 * np.array(luminosity_err)**2), 5)))
print("run_counter =", run_counter)

print("N_signal total =", reduce(lambda x, y: x + y, n_signal, 0))
print("N_background total =", reduce(lambda x, y: x + y, n_bckg, 0))
