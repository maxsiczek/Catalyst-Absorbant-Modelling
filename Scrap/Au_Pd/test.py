from ase.io import vasp
import os
from ase.build import fcc111, add_adsorbate # ASE's utilities to build the surface
from clusterx.parent_lattice import ParentLattice
from clusterx.structures_set import StructuresSet
from clusterx.visualization import juview
from clusterx.super_cell import SuperCell
from random import randint
from ase.io import read, write

# Build a 3-layer Al slab with vacancy in on-top configuration
pri2 = fcc111('Pd', size=(1, 1, 3))  # 3-atomic-layer Al slab
pri2.center(vacuum=10.0, axis=2)  # add vacuum along z-axis
platt2 = ParentLattice(pri2, site_symbols=[['Pd', 'Au'], ['Pd', 'Au'], ['Pd', 'Au']])
scell2 = SuperCell(platt2, [4, 4])
sset2 = StructuresSet(platt2)

num_files=2
Elist= []
Plist= []
for i in range(0,num_files):
    os.chdir(str({}).format(i))
    coords=vasp.read_vasp_out('OUTCAR')
    sset2.add_structure(coords)
    os.chdir('..')
    # print(atoms_object.get_cell())
    Elist.append(coords.get_potential_energy())
    Plist.append(coords.get_positions())
#print(Elist)
#print(str(len(Elist))+" energies")
#print(sset2)
#print(coords.get_positions(1))
#print(Plist)
#print(str(len(Plist))+" positions")
write('slab.xyz',sset2)
a=read('slab.xyz')
print(a)
