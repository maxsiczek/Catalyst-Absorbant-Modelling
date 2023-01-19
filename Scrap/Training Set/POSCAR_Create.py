#Creates POSCAR files from a specified number of random structures of 1 concentration of Au and Pd atoms


import os
from ase.build import fcc111, add_adsorbate  # ASE's utilities to build the surface
from clusterx.parent_lattice import ParentLattice
from clusterx.structures_set import StructuresSet
from clusterx.visualization import juview
from clusterx.super_cell import SuperCell
from random import randint
from ase.build import fcc111, add_adsorbate
from clusterx.super_cell import SuperCell
from clusterx.structures_set import StructuresSet
from ase.io import vasp
import numpy as np

nstruc=20

for i in range(0,nstruc):
    #Function to create folders
    def createFolder(directory):
        try:
            if not os.path.exists(directory):
                os.makedirs(directory)
        except OSError:
            print('Error: Creating directory. ' + directory)

    #Folder is created
    createFolder('./{}/'.format(i))

    pri2 = fcc111('Pd', size=(1, 1, 3))  # 3-atomic-layer PD slab
    pri2.center(vacuum=10.0, axis=2)  # add vacuum along z-axis

    # Build a 3-layer Pd and Au slab
    platt2 = ParentLattice(pri2, site_symbols=[['Pd', 'Au'], ['Pd', 'Au'], ['Pd', 'Au']])

    # Build the parent lattice and a 4x4 supercell
    scell2 = SuperCell(platt2, [4, 4])
    sset2 = StructuresSet(platt2)

    # concentration = {0:[randint(1,4*4)]} This is commented code to change concentration

    # Generate and add a random structure to the StructuresSet
    atoms_object = scell2.gen_random(24)
    sset2.add_structure(atoms_object)

    #Writes the structure to a POSCAR files in the a folder
    os.chdir(str(i))
    vasp.write_vasp('POSCAR', atoms_object,None,False,True)
    os.chdir('..')




