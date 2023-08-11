import subprocess as sub


energy_points = ["501", "503", "505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514", "517", "520", "525", "530"]
# energy_points = ["508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514"]
energy_points = ["509"]
aux1, aux2 = "\"", "\\"

for en in energy_points:
    command = f"root -l -q \"C:/work/Science/BINP/Kaon Mass Measure/up to date scripts for cutting/cutters/TrackEffCutter.cpp(\\{aux1 + en + aux2}\")\""
    res = sub.run(command, capture_output=True)
    output = res.stderr if res.stderr else res.stdout
    print(res.stdout.decode()[118:])
    print(res.stderr.decode())