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
sset6 = StructuresSet(platt2)



with open('structure.pkl', 'rb') as inp:
    strucs = pickle.load(inp)
    print(strucs)

sset6.add_structures(strucs)
print('first structure=',sset6.get_structure(0).get_sigmas())

from clusterx.calculators.emt import EMT2 # Load the EMT calculator from ASE
from clusterx.visualization import plot_property_vs_concentration

sset6.set_calculator(EMT2())
#temp=sset6.calculate_property("total_energy_emt") # Calculate energies with Effective Medium Theory calculator of ASE
temp=[-191.98633785, -192.3803963, -192.21625354, -191.81351726, -192.03387491, -191.69649921, -192.95811824, -192.46122795, -192.62040655, -191.66862616, -191.86058253, -191.84299077, -191.63256706, -191.70564561, -191.95900499, -192.29701507, -192.84809018, -192.31818476, -192.55573654, -191.37371426, -192.00352076, -193.14603609, -192.315369, -192.09209247, -191.61061585, -192.87835108, -191.78525661, -191.93109543, -191.95479766, -193.11085655, -192.46776992, -192.52801109, -192.34494702, -192.26066095, -192.22164227, -191.98759848, -191.86370915, -192.01021473, -192.34673888, -192.34437747, -190.73267162, -192.10860944, -192.16730105, -192.16601676, -192.94393568, -192.06158159, -191.86632121, -192.58051297, -191.81805869, -192.13172587, -191.98665716, -192.73453335, -191.23447735, -192.36955423, -191.35106556, -192.33892024, -191.92893451, -192.25113957, -192.47947183, -192.39863725, -191.99084895, -191.75831055, -192.69321787, -191.85313512, -191.97554422, -191.85180319, -192.6493575, -192.33859352, -191.94154621, -192.63414555, -191.1949463, -192.25953847, -191.98065361, -191.66115978, -191.58191491, -192.34619852, -191.34739675, -192.50454671, -192.27394753, -191.84511429, -192.02740716, -192.03305301, -192.01588955, -191.61178477, -191.84356671, -192.39518835, -192.01620164, -192.2660302, -192.49671185, -192.38635136, -192.98184883, -191.83210426, -192.01896905, -191.5939992, -191.29941827, -192.44401418, -191.17506294, -192.5299391, -191.87603351, -190.8848521]


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
from clusterx.monte_carlo import MonteCarloTrajectory

platt2.get_idx_subs()
mc1=MonteCarlo(cemodel1,scell2,{0:[46, 79]},no_of_swaps=1)
mc2=mc1.metropolis([8.6*10**-5,1000],10000,[46, 79, 79, 46, 46, 46, 79, 46, 79, 46, 79, 79, 46, 79, 79, 79, 46, 79, 46, 79, 46, 79, 46, 79, 46, 46, 46, 79, 79, 46, 79, 46, 46, 79, 79, 46, 79, 79, 46, 46, 46, 46, 79, 79, 79, 79, 46, 46])
nid=mc2.get_sampling_step_nos()
print('trajectory entries',nid)
print('Number of accepted=',len(mc2.get_model_total_energies()))


import matplotlib.pyplot as plt
plt.plot(mc2.get_model_total_energies()-mc2.get_model_total_energies()[0])
plt.xlabel("Number of Accepted Samples",fontsize=24)
plt.ylabel("Energy (eV)",fontsize=24)
# plt.title('Montecarlo Metropolis Sampling')
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.show()
# plt.hist(mc2.get_model_total_energies()-mc2.get_model_total_energies()[0])
# plt.xlabel("Energy (eV)",fontsize=12)
# plt.ylabel("Frequency",fontsize=12)
# plt.show()

# # vasp.write_vasp('POSCAR',mc2.get_structure(0))
# g0=np.zeros((len(nid),1))
# g1=np.zeros((len(nid),1))
# g2=np.zeros((len(nid),1))
# g3=np.zeros((len(nid),1))
#
# for k,i in enumerate(nid):
#     s=mc2.get_structure(i).get_atomic_numbers()
#     # unique, counts = np.unique(s, return_counts=True)
#     # print(dict(zip(unique, counts))) Tests concentration of the samples
#     t1=[s[2],s[5],s[14]]
#     t2=[s[8],s[17],s[20]]
#     t3=[s[26],s[29],s[38]]
#     t4=[s[35],s[47],s[44]]
#     for j in [t1,t2,t3,t4]:
#         if j.count(79) ==0:
#             g0[k]=np.sum(g0[0:k])+1)
#         if j.count(79) ==1:
#             g1[k]=np.sum(g1[0:k])+1)
#         if j.count(79) == 2:
#             g2[k]=np.sum(g2[0:k])+1)
#         if j.count(79) == 3:
#             g3[k]=np.sum(g3[0:k])+1)
# print('space')
# # print(g0/(4*len(nid)))
# # print(g1/(4*len(nid)))
# # print(g2/(4*len(nid)))
# # print(g3/(4*len(nid)))
# plt.plot(g0)
# plt.plot(g1)
# plt.plot(g2)
# plt.plot(g3)
# plt.show()

