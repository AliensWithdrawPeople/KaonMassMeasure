import subprocess as sub
from itertools import dropwhile
import re

energy_points = ["501", "503", "505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514"]
# energy_points = ["508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514"]
energy_points = ["514"]
aux1, aux2 = "\"", "\\"
pattern = r'luminosity = ([\d.]+)\r\nluminosity_err = ([\d.]+)\r\nentries = (\d+)'
        
luminosity = []
luminosity_err = []
n_events = []

for en in energy_points:
    command = f"root -l -q \"up to date scripts for cutting/cutter.cpp(\\{aux1 + en + aux2}\")\""
    res = sub.run(command, capture_output=True)
    output = res.stderr if res.stderr else res.stdout
    print(res.stdout.decode()[118:])
    
    match = re.search(pattern, res.stdout.decode())
    if match:
        lumi = float(match.group(1))
        lumi_err = float(match.group(2))
        events = int(match.group(3))
        n_events.append(events)
        luminosity.append(lumi)
        luminosity_err.append(lumi_err)

print("How many matches?", len(n_events)) 
print("Is it ok?", "Yes!" if (len(n_events) == len(energy_points)) else "No...")
print("n_events =", n_events)
print("luminosity =", luminosity)
print("luminosity_err =", luminosity_err)