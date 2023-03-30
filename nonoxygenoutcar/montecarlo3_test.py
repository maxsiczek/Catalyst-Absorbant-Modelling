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



with open('structure_no_oxygen.pkl', 'rb') as inp:
    strucs = pickle.load(inp)
    print(strucs)

sset6.add_structures(strucs)
print('first structure=',sset6.get_structure(0).get_sigmas())

from clusterx.calculators.emt import EMT2 # Load the EMT calculator from ASE
from clusterx.visualization import plot_property_vs_concentration

sset6.set_calculator(EMT2())
#temp=sset6.calculate_property("total_energy_emt") # Calculate energies with Effective Medium Theory calculator of ASE
temp=[-192.64947213, -192.22250837, -192.33364332, -192.25227456, -191.84170066, -192.1310434, -192.13449204, -192.05448528, -191.89050875, -192.04745238, -191.95556471, -191.94254732, -192.8753106, -192.99646496, -191.84211376, -191.60868848, -192.34065373, -191.95844178, -192.48282519, -192.25211357, -192.96175073, -192.01300452, -191.10478619, -191.61968718, -191.20984624, -192.45953042, -192.44510406, -191.52608294, -192.61563919, -192.85308486, -192.02095283, -191.9626645, -191.65174063, -192.37738382, -192.12413639, -191.49837309, -191.93507314, -192.45463796, -192.4559693, -192.57346049, -192.14214026, -192.25215371, -192.12256829, -192.75368256, -191.55719723, -191.72188046, -191.19747162, -191.25084179, -192.16154752, -191.88272261, -192.48627651, -192.2672662, -191.76589481, -192.47460761, -191.97112479, -192.2227161, -193.25304677, -192.36966571, -192.37165343, -191.63640973, -191.74317386, -192.4840845, -192.1523178, -191.26079156, -193.24923438, -191.85318116, -192.06647169, -191.61452741, -192.37045668, -191.39089627, -191.60422282, -191.42221292, -191.16978818, -192.4207627, -192.15469215, -192.24347073, -191.56172867, -192.25620512, -192.41322681, -192.03282623, -191.2106675, -191.81973175, -191.61886359, -191.96610437, -192.21548594, -192.87616701, -192.36033651, -191.98486517, -191.52151152, -192.09938557, -192.69932052, -191.47449131, -192.08331785, -192.44591596, -192.30435355, -192.30594017, -192.85653867, -192.63877354, -192.08050832, -192.21283578]


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
