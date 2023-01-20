from ase.io import vasp
from ase.constraints import FixAtoms
import os


for i in range(0,99):
    os.chdir(str(i))
    atoms_obj=vasp.read_vasp('POSCAR')
    c = FixAtoms(indices=[atom.index for atom in atoms_obj if atoms_obj.get_positions()[atom.index][2] < 11])
    atoms_obj.set_constraint(c) 
    vasp.write_vasp('POSCAR',atoms_obj,sort=True)
    os.chdir('..')
