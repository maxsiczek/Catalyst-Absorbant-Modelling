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


# list_of_atoms=[['Pd','Au'], ['Pd','Au'], ['Pd','Au'], ['X', 'O']]
# training_set_size=100
# concentration=[60,40]

class TsetNoAdsorbate:
    def __init__(self, list_of_atoms, training_set_size, concentration,directory_name):
        self.atoms = list_of_atoms
        self.size = training_set_size
        self.concentration = concentration
        self.directory=directory_name

    def sset_object(self):
        # Code that takes self.atoms and self.size to make the desired training set and exports a pickle file
        pri2 = fcc111(self.atoms[0], size=(1, 1, 3))  # 3-atomic-layer Pd slab
        pri2.center(vacuum=10.0, axis=2)  # add vacuum along z-axis
        # Build a 3-layer Au and Pd slab
        platt2 = ParentLattice(pri2, site_symbols=[self.atoms, self.atoms, self.atoms])
        # Build the parent lattice and a 4x4 supercell
        scell2 = SuperCell(platt2, [4, 4])
        sset = StructuresSet(platt2)

        nstruc = self.size

        for i in range(0, nstruc):
            atoms_object = scell2.gen_random(24)  # Makes random structure
            sset.add_structure(atoms_object)
        print(sset.get_structure(0).get_atomic_numbers())
        print(sset.get_nstr())
        print(sset.get_structure(0).get_atoms)
        print(sset.get_structure(0).get_supercell())

        strucs = sset.get_structures()
        def save_object(obj, filename):
            with open(filename, 'wb') as outp:  # Overwrites any existing file.
                pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)

        save_object(strucs, 'structure.pkl')

    def create_poscar_files(self):
        os.chdir(self.directory)
        pri2 = fcc111(self.atoms[0], size=(1, 1, 3))  # 3-atomic-layer Pd slab
        pri2.center(vacuum=10.0, axis=2)  # add vacuum along z-axis
        # Build a 3-layer Au and Pd slab
        platt2 = ParentLattice(pri2, site_symbols=[self.atoms, self.atoms, self.atoms])
        # Build the parent lattice and a 4x4 supercell
        scell2 = SuperCell(platt2, [4, 4])
        sset2 = StructuresSet(platt2)

        with open('structure.pkl', 'rb') as inp:
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
            atoms_object = scell2.gen_random(24)  # Makes random structure
            # concentration = {0:[randint(1,4*4)]}

            # Writes POSCAR file
            os.chdir(str(i))
            vasp.write_vasp('POSCAR', sset2.get_structure(i), None, False, True)
            os.chdir('../..')
        for i in range(0, self.size):
            files = ['INCAR', 'KPOINTS', 'POTCAR', 'relax.qs']
            for f in files:
                shutil.copy(f, '{}'.format(i))
        # Uses the training set found in the pickle file to create a POSCAR file for each structure in the training set

class TsetAdsorbate:
    def __init__(self, list_of_atoms, training_set_size, concentration,adsorbate,num_adsorbate,position):

        self.num_adsorbate=num_adsorbate
        self.size=training_set_size
        self.concentration=concentration

        pri2 = fcc111(list_of_atoms[0], size=(1, 1, 3))
        add_adsorbate(pri2, 'X', 1.7, position)  # on-top vacancy site
        pri2.center(vacuum=10.0, axis=2)  # add vacuum along z-axis

        print("{0:<19s}|{1:<19s}|{2:<19s}".format("Atom index", "Chemical symbol", "z coordinate"))
        for i, (symbol, z_coord) in enumerate(zip(pri2.get_chemical_symbols(), pri2.get_positions()[:, 2])):
            print("{0:<19d}|{1:<19s}|{2:<19.3f}".format(i, symbol, z_coord))
        platt2 = ParentLattice(pri2, site_symbols=[['Pd', 'Au'], ['Pd', 'Au'], ['Pd', 'Au'], ['X', adsorbate]])
        self.scell2 = SuperCell(platt2, [4, 4])
        self.scell2.get_sublattice_types(pretty_print=True)
        print(platt2.get_idx_subs())
        self.sset = StructuresSet(platt2)

    def sset_obj(self):
        sset1 = self.sset

        for i in range(self.size):
            sset1.add_structure(self.scell2.gen_random({0: [self.num_adsorbate], 1: [24]}))
        print(sset1.get_structure(0).get_atomic_numbers())
        print(sset1.get_nstr())
        print(sset1.get_structure(0).get_atoms)
        print(sset1.get_structure(0).get_supercell())

        strucs1 = sset1.get_structures()

        def save_object(obj, filename):
            with open(filename, 'wb') as outp:  # Overwrites any existing file.
                pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)

        save_object(strucs1, 'structure_oxygen.pkl')

    def create_poscar_files(self,potcar, ssetf):

        sset2 = self.sset

        with open(ssetf, 'rb') as inp:
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
            os.rename(potcar, "POTCAR")
            os.chdir('..')
    # create_poscar_files("POTCAR_oxygen",'structure_oxygen.pkl')
    def prepPoscar(self):
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
            replace_line("POSCAR", 0, 'Au  O Pd\n')
            replace_line("POSCAR", 5, ' Au  O   Pd\n')
            replace_line("POSCAR", 6, '  24   8  24\n')
            os.chdir('../..')