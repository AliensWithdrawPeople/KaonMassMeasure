import subprocess as sub
from itertools import dropwhile

energy_points = ["501", "503", "505", "508", "508.5", "509", "509.5", "510", "510.5", "511", "511.5", "514"]
aux1, aux2 = "\"", "\\"

for en in energy_points:
    command = f"root -l -q cutter.cpp(\\{aux1 + en + aux2}\")"
    res = sub.run(command, capture_output=True)
    output = res.stderr if res.stderr else res.stdout
    print(output)