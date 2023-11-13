import subprocess as sub
import re

energy_points = ["501", "503", "505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514", "517", "520", "525", "530"]

aux1, aux2 = "\"", "\\"

eff = []
mean_energies = []
mean_energies_err = []

output_pattern = r'Mean Energy = ([\d.]+) \+/- ([\d.]+); e_mc = ([\d.]+)'

def run_root(command: str):
    res = sub.run(command, capture_output=True)
    output = res.stderr if res.stderr else res.stdout    
    print(res.stdout.decode()[118:])
    print("ERROR!", res.stderr.decode())
    

for en in energy_points:
    command = f"root -l -q \"C:/work/Science/BINP/Kaon Mass Measure/up to date scripts for cutting/cutters/Xsec_EnergySpectrum_Cutter.cpp(\\{aux1 + en + aux2}\")\""
    run_root(command)
    
    command = f"root -l -q \"C:/work/Science/BINP/Kaon Mass Measure/up to date scripts for cutting/cutters/massCutter.cpp(\\{aux1 + en + aux2}\")\""
    run_root(command)

    
#     match = re.search(output_pattern, res.stdout.decode())
#     if match:
#         mean_energies.append(float(match.group(1)))
#         mean_energies_err.append(float(match.group(2)))
#         eff.append(float(match.group(3)))
    
# print("Means of Ks Energy Spectrums =", mean_energies)
# print("Mean errors of Ks Energy Spectrums =", mean_energies_err)
# print("e_MC =", eff)