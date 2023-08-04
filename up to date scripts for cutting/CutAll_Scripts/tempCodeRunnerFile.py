    match = re.search(output_pattern, res.stdout.decode())
#     if match:
#         mean_energies.append(float(match.group(1)))
#         mean_energies_err.append(float(match.group(2)))
#         eff.append(float(match.group(3)))
    
# print("Means of Ks Energy Spectrums =", mean_energies)
# print("Mean errors of Ks Energy Spectrums =", mean_energies_err)
# print("e_MC =",