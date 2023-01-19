from ase.io import vasp
import os
num_files=21
Plist= []
for i in range(0,num_files):
    os.chdir(str({}).format(i))
    coords=vasp.read_vasp_out('OUTCAR')
    os.chdir('..')
    # print(atoms_object.get_cell())
    Plist.append(coords.get_positions())
print(Plist)
print(str(len(Plist))+" positions")
