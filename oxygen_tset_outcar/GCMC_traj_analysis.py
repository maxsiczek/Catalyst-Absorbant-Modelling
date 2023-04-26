import ase
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
import clusterx
import random
import matplotlib.pyplot as plt
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





gl0 = []
gl1 = []
gl2 = []
gl3 = []
oxygencountlist= []

gb0 = []
gb1 = []
gb2 = []
gb3 = []
totalOxygenList = []

with open('GCMC_improved_trajectory.pkl', 'rb') as f:
    trajectory = pickle.load(f)
# print('test',trajectory[0].get_atomic_numbers())
atoms_list=[]
for obj in trajectory:
    atoms_list.append(obj.get_atoms())

ase.visualize.view(atoms_list)

for i in trajectory:
    g0 = 0
    g1 = 0
    g2 = 0
    g3 = 0
    s=i.get_atomic_numbers()
    # unique, counts = np.unique(s, return_counts=True)
    # print(dict(zip(unique, counts))) Tests concentration of the samples
    oxygencount=s.tolist().count(8)
    t1=[s[2],s[18],s[6]]
    t2=[s[26],s[30],s[42]]
    t3=[s[50],s[38],s[54]]
    t4=[s[46],s[58],s[62]]
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
    oxygencountlist.append(oxygencount)

    gl0.append(sum(gb0))
    gl1.append(sum(gb1))
    gl2.append(sum(gb2))
    gl3.append(sum(gb3))
    totalOxygenList.append(sum(oxygencountlist))



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
for i in range(1,len(trajectory)+1):
    rang.append(i*4)

totalVacancylist = [num for num in range(16, (len(trajectory)+1)*16, 16)]
oxygenPercentList = [num / totalVacancylist[i] for i, num in enumerate(totalOxygenList)]
gl0 = [int(b) / int(m) for b, m in zip(gl0, rang)]
gl1 = [int(b) / int(m) for b, m in zip(gl1, rang)]
gl2 = [int(b) / int(m) for b, m in zip(gl2, rang)]
gl3 = [int(b) / int(m) for b, m in zip(gl3, rang)]

print(gl0)
print(gl1)
print(gl2)
print(gl3)
print(oxygenPercentList)





print('space')
plt.plot(gl0,linewidth=3,color='red')
plt.plot(gl1,linewidth=3,color='blue')
plt.plot(gl2,linewidth=3,color='green')
plt.plot(gl3,linewidth=3,color='purple')
plt.plot(oxygenPercentList,linewidth=3,color='orange')
legend_drawn_flag = True
g0="Pd3"
g1="Au1Pd2"
g2="Au2Pd1"
g3="Au3"
oxygenPercentList="O Percent"

SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
plt.legend([g0.translate(SUB), g1.translate(SUB),g2.translate(SUB),g3.translate(SUB),oxygenPercentList], loc='upper right', frameon=legend_drawn_flag,fontsize=24)
plt.ylabel('Active Site Distribution (%)',fontsize=24)
plt.xlabel("Number of Accepted Samples",fontsize=24)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
# plt.title('Distribution of Active sites per Accepted Sample')
plt.show()
