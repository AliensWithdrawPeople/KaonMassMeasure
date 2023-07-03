import subprocess as sub
import uproot as up
import numpy as np

energy_points = ["501", "503", "505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514"]
aux1, aux2 = "\"", "\\"
eff = []

# for en in energy_points:
#     command = f"root -l -q \"up to date scripts for cutting/cutter_MC.cpp(\\{aux1 + en + aux2}\")\""
#     res = sub.run(command, capture_output=True)
#     output = res.stderr if res.stderr else res.stdout    
#     print(str(res.stdout))
#     print("ERROR!", str(res.stderr))


for en in energy_points:
    print(en)
    with up.open(f"C:/work/Science/BINP/Kaon Mass Measure/tr_ph/PhiXSection/MC/MC_Kn{en}.root:Kn") as ksTree: # type: ignore
        tree = ksTree.arrays(['emeas'], library="np") # type: ignore
    N_event1 = len(tree['emeas'])
    
    with up.open(f"C:/work/Science/BINP/Kaon Mass Measure/tr_ph/mcgpj/tr_ph v9/EnergySmearing/MCGPJ_kskl{en}_Merged.root:tr_ph_merged") as ksTree: # type: ignore
        tree = ksTree.arrays(['emeas'], library="np") # type: ignore
    N_event2 = len(tree['emeas'])
    eff.append(N_event1 / N_event2)    
print(list(np.round(eff, 5)))