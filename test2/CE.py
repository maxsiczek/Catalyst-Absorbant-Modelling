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
print(sset6.get_structure(0).get_sigmas())

from clusterx.calculators.emt import EMT2 # Load the EMT calculator from ASE
from clusterx.visualization import plot_property_vs_concentration

sset6.set_calculator(EMT2())
temp=sset6.calculate_property("total_energy_emt") # Calculate energies with Effective Medium Theory calculator of ASE
#temp=[-192.08582071, -191.87351734, -191.60018779, -191.39685004, -192.24290787, -191.78014638, -191.70745162, -192.97458495, -192.86267872, -191.37195796, -192.3143272, -192.14009737, -192.82585879, -192.05420399, -192.69999995, -192.47746079, -192.39000526, -192.26897121, -192.3791746, -192.43099895, -191.77825841, -191.61704212, -191.47478811, -192.03389507, -192.07830785, -192.69342232, -192.24191182, -192.54062519, -191.91839647, -192.51096022, -192.69325283, -192.17466602, -192.59825086, -191.72542479, -191.58574526, -192.14693661, -191.82264476, -192.03231015, -192.22406236, -192.13476996, -192.37122734, -192.14574102, -192.44934157, -192.40028879, -192.44270748, -192.91648491, -191.69067766, -191.88129597, -192.55731336, -191.98958383, -192.47554187, -192.06452466, -192.08804602, -192.6761159, -192.84187982, -192.76817289, -192.77311024, -191.68482559, -191.26662593, -191.75285029]


sset6.set_property_values(property_name='total_energy_emt', property_vals=temp)
print(temp)
r=3.5
from clusterx.clusters.clusters_pool import ClustersPool
cpool = ClustersPool(platt2, npoints=[0,1,2,3,4], radii=[0,0,r,r,r])
print(len(cpool)," clusters were generated.")
from clusterx.model import ModelBuilder



mb = ModelBuilder(selector_type="linreg",selector_opts={'clusters_sets':'size'},estimator_type="skl_LinearRegression",estimator_opts={"fit_intercept":False})
cemodel1 = mb.build(sset6, cpool, "total_energy_emt") #Build CE model using the training data set
cpool_opt1 = mb.get_opt_cpool()


cemodel1.report_errors(sset6)
cpool_opt1.display_info(ecis=cemodel1.get_ecis())

from clusterx.visualization import plot_optimization_vs_number_of_clusters
from clusterx.visualization import plot_predictions_vs_target
plot_optimization_vs_number_of_clusters(mb.get_selector(),scale=0.7)
plot_predictions_vs_target(sset6,cemodel1,"total_energy_emt",scale=0.7)