import csv
import math
import sys
from tqdm import tqdm
import numpy as np

ELEMENTARY_CHARGE = 1.602176634e-12 #erg
V0 = 2.19e8 # cm/s
a0 = 0.529e-8 # сm
IH = 13.6 # eV
da = 1.6605e-24 #g
piko2centi = 1e-10

species = input("species: ")
max_energy = float(input("max energy, eV: "))
energy_step = float(input("energy step, eV: "))
source_mass = float(input("source mass, dalton: "))
target_mass = float(input("target mass, dalton: "))
gk_radius = float(input("gk radius, pm: "))
polarization_alpha_per_a03 = float(input("polarizability a.u.: "))
I = float(input("ionization threshold, eV: "))
cs_polarization_elastic = []
cs_charge_exchange = []
print("calculating")
for i in tqdm(np.geomspace(1,int(max_energy/energy_step),100)):
    energy = energy_step*i
    cs_pol_qm = 2*math.sqrt(2)*math.pi*a0*a0*math.sqrt(polarization_alpha_per_a03*IH/energy)
    cs_pol_hs = 0#math.pi*(gk_radius*piko2centi)**2
    cs_polarization_elastic.append((energy, 1e-4*math.sqrt(cs_pol_qm**2+cs_pol_hs**2)))

    v_rel_ion = math.sqrt(2*energy*ELEMENTARY_CHARGE/(da*source_mass)) #cm/s
    cs_charge_exchange.append((energy, 1e-4*math.pi*a0*a0*IH/I*math.pow(math.log(100*V0/v_rel_ion * math.sqrt(I/IH)),2))) #m^2


print("writing")
with open(f"{species}_Raizer.txt", 'w') as file:
    file.write("ION_BACKSCATTER\n")
    file.write(f"{species}^+ / {species} -> {species}^+\n")
    file.write(f"{source_mass/target_mass} {target_mass}\n")
    file.write("---------\n")
    file.write(f"{0.0}	{0.0}\n")
    for energy,cs in tqdm(cs_charge_exchange):
        file.write(f"{energy}	{cs}\n")
    file.write("---------\n")

    file.write("\n\n\n")

    file.write("ION_ISOTROPIC\n")
    file.write(f"{species}^+ / {species} -> {species}^+\n")
    file.write(f"{source_mass/target_mass} {target_mass}\n")
    file.write("---------\n")
    file.write(f"{0.0}	{0.0}\n")
    for energy,cs in tqdm(cs_polarization_elastic):
        file.write(f"{energy}	{cs}\n")
    file.write("---------\n")