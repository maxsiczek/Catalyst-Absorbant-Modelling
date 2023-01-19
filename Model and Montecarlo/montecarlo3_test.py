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
mc2=mc1.metropolis([8.6*10**-5,1000],20000,[46, 79, 79, 46, 46, 46, 79, 46, 79, 46, 79, 79, 46, 79, 79, 79, 46, 79, 46, 79, 46, 79, 46, 79, 46, 46, 46, 79, 79, 46, 79, 46, 46, 79, 79, 46, 79, 79, 46, 46, 46, 46, 79, 79, 79, 79, 46, 46])
nid=mc2.get_sampling_step_nos()
print('trajectory entries',nid)
print('Number of accepted=',len(mc2.get_model_total_energies()))


import matplotlib.pyplot as plt
# plt.plot(mc2.get_model_total_energies()-mc2.get_model_total_energies()[0])
# plt.xlabel("Number of Sampling Steps")
# plt.ylabel("Energy (eV)")
# plt.show()
# plt.hist(mc2.get_model_total_energies()-mc2.get_model_total_energies()[0])
# plt.xlabel("Energy (eV)")
# plt.ylabel("Frequency")
# plt.show()

vasp.write_vasp('POSCAR',mc2.get_structure_at_step(1000))

gl0 = []
gl1 = []
gl2 = []
gl3 = []

gb0 = []
gb1 = []
gb2 = []
gb3 = []



for i in nid:
    g0 = 0
    g1 = 0
    g2 = 0
    g3 = 0
    s=mc2.get_structure_at_step(i).get_atomic_numbers()
    # unique, counts = np.unique(s, return_counts=True)
    # print(dict(zip(unique, counts))) Tests concentration of the samples
    t1=[s[2],s[5],s[14]]
    t2=[s[8],s[17],s[20]]
    t3=[s[26],s[29],s[38]]
    t4=[s[35],s[47],s[44]]
    for j in [t1, t2, t3, t4]:
        if j.count(79) == 0:
            g0 = g0 + 1
        if j.count(79) == 1:
            g1 = g1 + 1
        if j.count(79) == 2:
            g2 = g2 + 1
        if j.count(79) == 3:
            g3 = g3 + 1
    gb0.append(g0)
    gb1.append(g1)
    gb2.append(g2)
    gb3.append(g3)

    gl0.append(sum(gb0))
    gl1.append(sum(gb1))
    gl2.append(sum(gb2))
    gl3.append(sum(gb3))



print(gb0)
print(gb1)
print(gb2)
print(gb3)

print('space')

print(gl0)
print(gl1)
print(gl2)
print(gl3)


rang=[]
for i in range(1,len(nid)+1):
    rang.append(i*4)


gl0 = [int(b) / int(m) for b, m in zip(gl0, rang)]
gl1 = [int(b) / int(m) for b, m in zip(gl1, rang)]
gl2 = [int(b) / int(m) for b, m in zip(gl2, rang)]
gl3 = [int(b) / int(m) for b, m in zip(gl3, rang)]

print(gl0)
print(gl1)
print(gl2)
print(gl3)





print('space')
plt.plot(gl0)
plt.plot(gl1)
plt.plot(gl2)
plt.plot(gl3)
legend_drawn_flag = True
g0="Pd3"
g1="Au1Pd2"
g2="Au2Pd1"
g3="Au3"

SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
plt.legend([g0.translate(SUB), g1.translate(SUB),g2.translate(SUB),g3.translate(SUB)], loc='upper right', frameon=legend_drawn_flag,fontsize=24)
plt.ylabel('Active Site Distribution (%)',fontsize=24)
plt.xlabel("Number of Accepted Samples",fontsize=24)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
# plt.title('Distribution of Active sites per Accepted Sample')
plt.show()


cutoff=18500

del gl0[:cutoff]
del gl1[:cutoff]
del gl2[:cutoff]
del gl3[:cutoff]

gm0=np.mean(gl0)
gm1=np.mean(gl1)
gm2=np.mean(gl2)
gm3=np.mean(gl3)

print(gm0)
print(gm1)
print(gm2)
print(gm3)



# plt.plot(gl0)
# plt.plot(gl1)
# plt.plot(gl2)
# plt.plot(gl3)
# plt.xlabel('Number of Structures')
# plt.ylabel('Active Site Distribution (%)')
#
# g0="Pd3"
# g1="Au1Pd2"
# g2="Au2Pd1"
# g3="Au3"
#
# SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
# legend_drawn_flag = True
# plt.legend([g0.translate(SUP), g1.translate(SUP),g2.translate(SUP),g3.translate(SUP)], loc=0, frameon=legend_drawn_flag)
# plt.show()
