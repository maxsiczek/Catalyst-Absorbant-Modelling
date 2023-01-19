#Extracts the energies from OUTCAR files in the same folder and makes them into a list. Number of OUTCAR files must be specified.


from ase.io import vasp
import os
num_files=21
Elist= []
for i in range(0,num_files):
    os.chdir(str({}).format(i))
    coords=vasp.read_vasp_out('OUTCAR')
    os.chdir('..')
    # print(atoms_object.get_cell())
    Elist.append(coords.get_potential_energy())
print(Elist)
print(str(len(Elist))+" energies")
