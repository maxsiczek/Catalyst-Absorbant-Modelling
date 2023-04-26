import ase.visualize
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
from clusterx import structure
from clusterx.calculators.emt import EMT2 # Load the EMT calculator from ASE
from clusterx.visualization import plot_property_vs_concentration
from ase.io import vasp
import os
import pickle
import clusterx
import random
from ase import *

pri2 = fcc111('Pd', size=(1, 1, 3))
add_adsorbate(pri2, 'X', 1.7, 'fcc')  # on-top vacancy site
pri2.center(vacuum=10.0, axis=2)  # add vacuum along z-axis

print("{0:<19s}|{1:<19s}|{2:<19s}".format("Atom index", "Chemical symbol", "z coordinate"))
for i, (symbol, z_coord) in enumerate(zip(pri2.get_chemical_symbols(), pri2.get_positions()[:, 2])):
    print("{0:<19d}|{1:<19s}|{2:<19.3f}".format(i, symbol, z_coord))
platt2 = ParentLattice(pri2, site_symbols=[['Pd', 'Au'], ['Pd', 'Au'], ['Pd', 'Au'], ['X', 'O']])
scell2 = SuperCell(platt2, [4, 4])
scell2.get_sublattice_types(pretty_print=True)
sset6 = StructuresSet(platt2)




with open('structure_oxygen.pkl', 'rb') as inp:
    strucs = pickle.load(inp)
    print(strucs)

sset6.add_structures(strucs)
print('first structure=',sset6.get_structure(0).get_atomic_numbers())

from clusterx.calculators.emt import EMT2 # Load the EMT calculator from ASE
from clusterx.visualization import plot_property_vs_concentration

sset6.set_calculator(EMT2())
#temp=sset6.calculate_property("total_energy_emt") # Calculate energies with Effective Medium Theory calculator of ASE
temp=[-233.75019001, -230.39949145, -230.41642605, -232.92550792, -232.25474271, -232.06625406, -231.03230612, -231.01794057, -232.53805623, -232.0541816, -230.56729887, -231.83706403, -231.64181822, -232.64229974, -231.1761779, -233.37888234, -232.29175595, -231.25529331, -231.14917998, -233.27754977, -230.51711417, -232.33493018, -230.56351581, -232.71507182, -232.21259049, -230.11796393, -230.52814574, -232.9856879, -232.2989704, -232.79855691, -229.58941518, -233.89702516, -231.90692845, -231.45051152, -234.68676429, -232.33541047, -232.32402045, -231.58213114, -233.91299888, -232.19509312, -232.05603928, -230.37099549, -232.88551241, -231.7170679, -230.4244832, -233.53623428, -233.60564537, -231.37999226, -230.42961338, -230.37825858, -232.15860339, -230.8508302, -231.11564108, -231.40353516, -230.67177557, -230.90346862, -233.29469413, -230.92520531, -233.31731473, -232.8347633, -230.34773993, -229.08039305, -230.31953305, -230.64803972, -232.59785396, -231.84562338, -233.80359322, -232.11151289, -231.77306496, -230.58876462, -228.45575599, -235.12152543, -235.0470921, -233.6964058, -230.49305883, -232.32509228, -230.07314787, -231.72602331, -228.67022594, -230.58062675, -232.05714829, -233.11207428, -231.78318343, -232.81874704, -231.69051811, -231.04156462, -232.13916759, -232.26584581, -232.1842805, -233.6457825, -233.16723773, -232.74395021, -232.81595047, -232.41401241, -231.75681963, -232.08839623, -229.31978215, -230.24534805, -230.86939434, -232.939754]


sset6.set_property_values(property_name='total_energy_emt', property_vals=temp)
print('energies=',temp)
r=3.5
from clusterx.clusters.clusters_pool import ClustersPool
cpool = ClustersPool(platt2, npoints=[0,1,2,3,4], radii=[0,0,r,r,r])
print(len(cpool)," clusters were generated.")
from clusterx.model import ModelBuilder



mb = ModelBuilder(selector_type="linreg",selector_opts={'clusters_sets':'size'},estimator_type="skl_LinearRegression",estimator_opts={"fit_intercept":False})
cemodel1 = mb.build(sset6, cpool, "total_energy_emt") #Build CE model using the training data set
cpool_opt1 = mb.get_opt_cpool()


# cemodel1.report_errors(sset6)
# cpool_opt1.display_info(ecis=cemodel1.get_ecis())
from clusterx.visualization import plot_optimization_vs_number_of_clusters
from clusterx.visualization import plot_predictions_vs_target
# plot_optimization_vs_number_of_clusters(mb.get_selector(),scale=0.7)
# plot_predictions_vs_target(sset6,cemodel1,"total_energy_emt",scale=0.7)

#Montecarlo
from clusterx.monte_carlo import MonteCarlo

N_total = 5000
init_struc = [79, 79, 79,  0, 79, 46, 46,  8 ,79, 46, 46,  0, 46, 79 ,46, 8, 46 ,79, 46 , 0,46, 46, 46 , 8,46, 79 ,79 , 0 ,79, 46, 79 , 8,46 ,79 ,79 , 0 ,79 ,79, 46 , 8,79 ,46, 46 , 0 ,79 ,79, 79,  8, 79 ,46 ,79 , 0 ,46 ,46 ,46,  8 ,79 ,46, 79,  0,46 ,46 ,79,  8]
GCMC_struc = init_struc
trajectory = []

# total num of vacancies
N_tot = 16
beta = (8.6 * 10 ** (-5) * 1000)
chem_pot = 10**0/ beta


for i in range(0, N_total):
    # choose with equal probability
    # 1, 2, or 3
    random_step = random.choice([1,2,3])


    # 1 swap MC (regular MC) move
    if random_step == 1:
        mc1 = MonteCarlo(cemodel1, scell2, {0: [46, 79], 1 : [0, 8] }, no_of_swaps=1,sublattice_indices=[1])
        mc2 = mc1.metropolis([8.6 * 10 ** -5, 1000], 1, GCMC_struc)
        nentries = len(mc2.get_model_total_energies())
        GCMC_struc=mc2.get_structure_at_step(nentries).get_atomic_numbers().tolist()
        print('Entries', nentries)
        if nentries == 2:
            trajectory.append(structure.Structure(scell2, decoration=GCMC_struc))
        else:
            GCMC_struc=GCMC_struc


    # 2 O addition GCMC move
    # count the num of O atoms in the structure set that equal to No
    if random_step == 2 and GCMC_struc.count(8) < N_tot:
        nO=GCMC_struc.count(8)
        # Randomly select one of the indices of the 0's in the list
        random_vacancy = random.choice([i for i, x in enumerate(GCMC_struc) if x == 0])

        # Create a new list that is a copy of the original list with the selected 0 changed to 8
        new_struc = [8 if i == random_vacancy else x for i, x in enumerate(GCMC_struc)]



        # randomly choose on of these
        # make a new strcutre with the O atom added^


        # calculate Energy of this new struc = E2
        # enregy of the original struc = E1
        print(GCMC_struc)
        print(new_struc)
        E1 = cemodel1.predict(structure.Structure(scell2, decoration=GCMC_struc))
        E2 = cemodel1.predict(structure.Structure(scell2, decoration=new_struc))

         # calc. acceptance ratio

        acc_ratio = ((N_tot - nO) / (nO + 1)) * np.exp(-beta * (E2 - E1) + beta * chem_pot)

    # draw a random num between 0 and 1 = rand_num_1
        rand_num_1 = random.random()

        if rand_num_1 < acc_ratio:
            GCMC_struc = new_struc
            print('addtion accepted')
            trajectory.append(structure.Structure(scell2, decoration=GCMC_struc))
        else:
            GCMC_struc = GCMC_struc
            print('addtion not accepted')


    if random_step == 3 and GCMC_struc.count(8) > 0 :
        nO=GCMC_struc.count(8)
        # Randomly select one of the indices of the 8's in the list
        random_vacancy = random.choice([i for i, x in enumerate(GCMC_struc) if x == 8])

        # Create a new list that is a copy of the original list with the selected 8 changed to 0
        new_struc = [0 if i == random_vacancy else x for i, x in enumerate(GCMC_struc)]


        # randomly choose on of these
        # make a new strucutre with O deleted


        # calculate Energy of this new struc = E2
        # enregy of the original struc = E1
        print(GCMC_struc)
        print(new_struc)
        E1 = cemodel1.predict(structure.Structure(scell2, decoration=GCMC_struc))
        E2 = cemodel1.predict(structure.Structure(scell2, decoration=new_struc))

         # calc. acceptance ratio

        acc_ratio = (nO / (N_tot - (nO - 1))) * np.exp(-beta * (E2 - E1) - beta * chem_pot)

    # draw a random num between 0 and 1 = rand_num_1
        rand_num_1 = random.random()

        if rand_num_1 < acc_ratio:
            GCMC_struc = new_struc
            print('deletion accepted')
            trajectory.append(structure.Structure(scell2, decoration=GCMC_struc))
        else:
            GCMC_struc = GCMC_struc
            print('deletion not accepted')


print(trajectory)
print(len(trajectory))
with open('GCMC_improved_trajectory.pkl', 'wb') as f:
    pickle.dump(trajectory, f)

