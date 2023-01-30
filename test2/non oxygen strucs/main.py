# import CreateTrainingSet
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
import shutil, os

def create_poscar_files(potcar,sset):
    pri2 = fcc111('Pd', size=(1, 1, 3))
    pri2.center(vacuum=10.0, axis=2)  # add vacuum along z-axis

    print("{0:<19s}|{1:<19s}|{2:<19s}".format("Atom index", "Chemical symbol", "z coordinate"))
    for i, (symbol, z_coord) in enumerate(zip(pri2.get_chemical_symbols(), pri2.get_positions()[:, 2])):
        print("{0:<19d}|{1:<19s}|{2:<19.3f}".format(i, symbol, z_coord))
    platt2 = ParentLattice(pri2, site_symbols=[['Pd', 'Au'], ['Pd', 'Au'], ['Pd', 'Au']])
    scell2 = SuperCell(platt2, [4, 4])
    scell2.get_sublattice_types(pretty_print=True)
    sset2 = StructuresSet(platt2)

    with open(sset, 'rb') as inp:
        strucs = pickle.load(inp)
        print(strucs)

    sset2.add_structures(strucs)
    nstruc = sset2.get_nstr()
    print(nstruc)
    print(sset2.get_structure(0).get_sigmas())

    # Function to create folders
    def createFolder(directory):
        try:
            if not os.path.exists(directory):
                os.makedirs(directory)
        except OSError:
            print('Error: Creating directory. ' + directory)

    for i in range(0, nstruc):
        createFolder('./{}/'.format(i))
        os.chdir(str(i))
        vasp.write_vasp('POSCAR', sset2.get_structure(i), sort=True)
        os.chdir('..')
    for i in range(0, nstruc):
        createFolder('Viewsset')
        os.chdir('Viewsset')
        vasp.write_vasp('POSCAR_{}'.format(i), sset2.get_structure(i))
        os.chdir('..')
    for i in range(0, nstruc):
        files = ['INCAR', 'KPOINTS', potcar, 'relax.qs']
        for f in files:
            shutil.copy(f, '{}'.format(i))
        os.chdir(str(i))
        os.rename(potcar,"POTCAR")
        os.chdir('..')
create_poscar_files('POTCAR',"structure_no_oxygen.pkl")