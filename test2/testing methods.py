import CreateTrainingSet
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
# module.CreateTrainingSet(['Pd','Au'],100,None).create_structure_set_object()
# module.CreateTrainingSet(['Pd','Au'],100,None).create_poscar_files()


def sset_with_oxygen():
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
    sset = StructuresSet(platt2)

    nstruc = 100

    for i in range(nstruc):
        sset.add_structure(scell2.gen_random({0:[8],1:[24]}))
    print(sset.get_structure(0).get_atomic_numbers())
    print(sset.get_nstr())
    print(sset.get_structure(0).get_atoms)
    print(sset.get_structure(0).get_supercell())

    strucs = sset.get_structures()

    def save_object(obj, filename):
        with open(filename, 'wb') as outp:  # Overwrites any existing file.
            pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)

    save_object(strucs, 'structure_oxygen.pkl')

# sset_with_oxygen()
# CreateTrainingSet.CreateTrainingSet(['Au','Pd'],100,None).create_poscar_files()


def create_poscar_files(potcar,sset):
    pri2 = fcc111('Pd', size=(1, 1, 3))
    add_adsorbate(pri2, 'X', 1.7, 'fcc')  # on-top vacancy site
    pri2.center(vacuum=10.0, axis=2)  # add vacuum along z-axis

    print("{0:<19s}|{1:<19s}|{2:<19s}".format("Atom index", "Chemical symbol", "z coordinate"))
    for i, (symbol, z_coord) in enumerate(zip(pri2.get_chemical_symbols(), pri2.get_positions()[:, 2])):
        print("{0:<19d}|{1:<19s}|{2:<19.3f}".format(i, symbol, z_coord))
    platt2 = ParentLattice(pri2, site_symbols=[['Pd', 'Au'], ['Pd', 'Au'], ['Pd', 'Au'], ['X', 'O']])
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
# create_poscar_files("POTCAR_oxygen",'structure_oxygen.pkl')

def prepPoscar():
    def replace_line(file_name, line_num, text):
        lines = open(file_name, 'r').readlines()
        lines[line_num] = text
        out = open(file_name, 'w')
        out.writelines(lines)
        out.close()

    for i in range(0, 100):
        os.chdir(str(i))
        with open(r"POSCAR", 'r+') as fp:
            # read an store all lines into list
            lines = fp.readlines()
            # move file pointer to the beginning of a file
            fp.seek(0)
            # truncate the file
            fp.truncate()
            fp.writelines(lines[:-8])
            # start writing lines
            # iterate line and line number
        replace_line("POSCAR",0,'Au  O Pd\n')
        replace_line("POSCAR", 5, ' Au  O   Pd\n')
        replace_line("POSCAR", 6, '  24   8  24\n')
        os.chdir('..')

# prepPoscar()