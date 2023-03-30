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
print(sset6.get_structure(0).get_sigmas())

from clusterx.calculators.emt import EMT2 # Load the EMT calculator from ASE
from clusterx.visualization import plot_property_vs_concentration

sset6.set_calculator(EMT2())
#temp=sset6.calculate_property("total_energy_emt") # Calculate energies with Effective Medium Theory calculator of ASE
temp=[-192.64947213, -192.22250837, -192.33364332, -192.25227456, -191.84170066, -192.1310434, -192.13449204, -192.05448528, -191.89050875, -192.04745238, -191.95556471, -191.94254732, -192.8753106, -192.99646496, -191.84211376, -191.60868848, -192.34065373, -191.95844178, -192.48282519, -192.25211357, -192.96175073, -192.01300452, -191.10478619, -191.61968718, -191.20984624, -192.45953042, -192.44510406, -191.52608294, -192.61563919, -192.85308486, -192.02095283, -191.9626645, -191.65174063, -192.37738382, -192.12413639, -191.49837309, -191.93507314, -192.45463796, -192.4559693, -192.57346049, -192.14214026, -192.25215371, -192.12256829, -192.75368256, -191.55719723, -191.72188046, -191.19747162, -191.25084179, -192.16154752, -191.88272261, -192.48627651, -192.2672662, -191.76589481, -192.47460761, -191.97112479, -192.2227161, -193.25304677, -192.36966571, -192.37165343, -191.63640973, -191.74317386, -192.4840845, -192.1523178, -191.26079156, -193.24923438, -191.85318116, -192.06647169, -191.61452741, -192.37045668, -191.39089627, -191.60422282, -191.42221292, -191.16978818, -192.4207627, -192.15469215, -192.24347073, -191.56172867, -192.25620512, -192.41322681, -192.03282623, -191.2106675, -191.81973175, -191.61886359, -191.96610437, -192.21548594, -192.87616701, -192.36033651, -191.98486517, -191.52151152, -192.09938557, -192.69932052, -191.47449131, -192.08331785, -192.44591596, -192.30435355, -192.30594017, -192.85653867, -192.63877354, -192.08050832, -192.21283578]


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

print(cemodel1.predict(sset6.get_structure(0)))#test



from clusterx.visualization import plot_optimization_vs_number_of_clusters
from clusterx.visualization import plot_predictions_vs_target
plot_optimization_vs_number_of_clusters(mb.get_selector(),scale=0.7)
plot_predictions_vs_target(sset6,cemodel1,"total_energy_emt",scale=0.7)


