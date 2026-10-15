from ase.io import read
import numpy as np

# 1. Wczytanie komórki bazowej z pliku CIF
unit_cell = read("Al2Cu.cif")

# 2. Powielenie komórki (tworzymy superkomórkę np. 6x6x6)
# Standardowa komórka Al2Cu ma 12 atomów. 12 * 4 * 4 * 5 = 2592 atomy
supercell = unit_cell.repeat((10, 10, 10))

# 3. Sortowanie: Twój kod w C++ wymaga, aby najpierw było Al, potem Cu
al_atoms = [atom for atom in supercell if atom.symbol == 'Al']
cu_atoms = [atom for atom in supercell if atom.symbol == 'Cu']

# 4. Zapis do prostego pliku tekstowego (tylko x, y, z)
with open("start_pos.dat", "w") as f:
    for atom in al_atoms + cu_atoms:
        f.write(f"{atom.position[0]:.6f} {atom.position[1]:.6f} {atom.position[2]:.6f}\n")

# Zdobądź parametry dla C++
N_calkowite = len(al_atoms) + len(cu_atoms)
Lx = supercell.cell.lengths()[0] # zakladamy ze pudelko jest szescienne L_x = L_y = L_z
Ly = supercell.cell.lengths()[1]
Lz = supercell.cell.lengths()[2]

print("GOTOWE!")
print(f"Wpisz w C++: params.N = {N_calkowite};")
print(f"Wpisz w C++: params.Lx = {Lx:.5f}; params.Ly = {Ly:.5f}; params.Lz = {Lz:.5f}; ")