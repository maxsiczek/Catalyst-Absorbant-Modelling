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
from clusterx.calculators.emt import EMT2 # Load the EMT calculator from ASE
from clusterx.visualization import plot_property_vs_concentration
from ase.io import vasp
import os
import pickle


pri2 = fcc111('Pd', size=(1, 1, 3))  # 3-atomic-layer Pd slab
pri2.center(vacuum=10.0, axis=2)  # add vacuum along z-axis
# Build a 3-layer Au and Pd slab
platt2 = ParentLattice(pri2, site_symbols=[['Pd', 'Au'], ['Pd', 'Au'], ['Pd', 'Au']])
# Build the parent lattice and a 4x4 supercell
scell2 = SuperCell(platt2, [4, 4])
sset2 = StructuresSet(platt2)



with open('structure.pkl', 'rb') as inp:
    strucs = pickle.load(inp)
    print(strucs)


sset2.add_structures(strucs)
nstruc=sset2.get_nstr()
print(nstruc)
print(sset2.get_structure(0).get_sigmas())




#Function to create folders
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

for i in range(0,nstruc):
    createFolder('./{}/'.format(i))
    atoms_object = scell2.gen_random(24) #Makes random structure
    # concentration = {0:[randint(1,4*4)]}

    #Writes POSCAR file
    os.chdir(str(i))
    vasp.write_vasp('POSCAR',sset2.get_structure(i) , None, False, True)
    os.chdir('..')



