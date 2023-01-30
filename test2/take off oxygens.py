import clusterx.structure
from ase.build import fcc111, add_adsorbate  # ASE's utilities to build the surface
from clusterx.parent_lattice import ParentLattice
from clusterx.structures_set import StructuresSet
from clusterx.visualization import juview
from clusterx.super_cell import SuperCell
from random import randint
from ase.io import read, write
from ase.build import fcc111, add_adsorbate
import numpy as np
from clusterx.super_cell import SuperCell
import numpy as np
from clusterx.structures_set import StructuresSet
from clusterx.calculators.emt import EMT2  # Load the EMT calculator from ASE
from clusterx.visualization import plot_property_vs_concentration
from ase.io import vasp
import os
import pickle
import shutil, os
from ase.io import vasp
from ase.constraints import FixAtoms




pri2 = fcc111('Pd', size=(1, 1, 3))
add_adsorbate(pri2, 'X', 1.7, 'fcc')  # on-top vacancy site
pri2.center(vacuum=10.0, axis=2)  # add vacuum along z-axis

print("{0:<19s}|{1:<19s}|{2:<19s}".format("Atom index", "Chemical symbol", "z coordinate"))
for i, (symbol, z_coord) in enumerate(zip(pri2.get_chemical_symbols(), pri2.get_positions()[:, 2])):
    print("{0:<19d}|{1:<19s}|{2:<19.3f}".format(i, symbol, z_coord))
platt2 = ParentLattice(pri2, site_symbols=[['Pd','Au'], ['Pd','Au'], ['Pd','Au'], ['X', 'O']])
scell2 = SuperCell(platt2, [4, 4])
scell2.get_sublattice_types(pretty_print=True)
print(platt2.get_idx_subs())
sset2 = StructuresSet(platt2)

with open('structure_oxygen.pkl', 'rb') as inp:
    strucs = pickle.load(inp)
    print(strucs)

sset2.add_structures(strucs)

pri1 = fcc111('Pd', size=(1, 1, 3))
pri1.center(vacuum=10.0, axis=2)  # add vacuum along z-axis

print("{0:<19s}|{1:<19s}|{2:<19s}".format("Atom index", "Chemical symbol", "z coordinate"))
for i, (symbol, z_coord) in enumerate(zip(pri1.get_chemical_symbols(), pri1.get_positions()[:, 2])):
    print("{0:<19d}|{1:<19s}|{2:<19.3f}".format(i, symbol, z_coord))
platt1 = ParentLattice(pri1, site_symbols=[['Pd','Au'], ['Pd','Au'], ['Pd','Au']])
scell1 = SuperCell(platt1, [4, 4])
scell1.get_sublattice_types(pretty_print=True)
print(platt2.get_idx_subs())
sset1 = StructuresSet(platt1)


for i in range(0,100):
    dec=sset2.get_structure(i).get_atomic_numbers()
    dec=np.delete(dec,[3,7,11,15,19,23,27,31,35,39,43,47,51,55,59,63])
    s=clusterx.structure.Structure(scell1,decoration=dec)
    print(s.get_atomic_numbers())
    sset1.add_structure(s)
def save_object(obj, filename):
    with open(filename, 'wb') as outp:  # Overwrites any existing file.
        pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)

struc=sset1.get_structures()
save_object(struc, 'non oxygen strucs/structure_no_oxygen.pkl')






