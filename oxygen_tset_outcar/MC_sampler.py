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
from clusterx.monte_carlo import MonteCarlo
from clusterx.monte_carlo import MonteCarloTrajectory

class MCsim:
    def __init__(self):
        pri2 = fcc111('Pd', size=(1, 1, 3))
        add_adsorbate(pri2, 'X', 1.7, 'fcc')  # on-top vacancy site
        pri2.center(vacuum=10.0, axis=2)  # add vacuum along z-axis

        print("{0:<19s}|{1:<19s}|{2:<19s}".format("Atom index", "Chemical symbol", "z coordinate"))
        for i, (symbol, z_coord) in enumerate(zip(pri2.get_chemical_symbols(), pri2.get_positions()[:, 2])):
            print("{0:<19d}|{1:<19s}|{2:<19.3f}".format(i, symbol, z_coord))
        platt2 = ParentLattice(pri2, site_symbols=[['Pd', 'Au'], ['Pd', 'Au'], ['Pd', 'Au'], ['X', 'O']])
        self.scell2 = SuperCell(platt2, [4, 4])
        self.scell2.get_sublattice_types(pretty_print=True)
        self.sset6 = StructuresSet(platt2)

        with open('structure_oxygen.pkl', 'rb') as inp:
            strucs = pickle.load(inp)
            print(strucs)

        self.sset6.add_structures(strucs)
        print('first structure=', self.sset6.get_structure(0).get_atomic_numbers())

        from clusterx.calculators.emt import EMT2  # Load the EMT calculator from ASE
        from clusterx.visualization import plot_property_vs_concentration

        self.sset6.set_calculator(EMT2())
        # temp=sset6.calculate_property("total_energy_emt") # Calculate energies with Effective Medium Theory calculator of ASE
        temp = [-233.75019001, -230.39949145, -230.41642605, -232.92550792, -232.25474271, -232.06625406, -231.03230612,
                -231.01794057, -232.53805623, -232.0541816, -230.56729887, -231.83706403, -231.64181822, -232.64229974,
                -231.1761779, -233.37888234, -232.29175595, -231.25529331, -231.14917998, -233.27754977, -230.51711417,
                -232.33493018, -230.56351581, -232.71507182, -232.21259049, -230.11796393, -230.52814574, -232.9856879,
                -232.2989704, -232.79855691, -229.58941518, -233.89702516, -231.90692845, -231.45051152, -234.68676429,
                -232.33541047, -232.32402045, -231.58213114, -233.91299888, -232.19509312, -232.05603928, -230.37099549,
                -232.88551241, -231.7170679, -230.4244832, -233.53623428, -233.60564537, -231.37999226, -230.42961338,
                -230.37825858, -232.15860339, -230.8508302, -231.11564108, -231.40353516, -230.67177557, -230.90346862,
                -233.29469413, -230.92520531, -233.31731473, -232.8347633, -230.34773993, -229.08039305, -230.31953305,
                -230.64803972, -232.59785396, -231.84562338, -233.80359322, -232.11151289, -231.77306496, -230.58876462,
                -228.45575599, -235.12152543, -235.0470921, -233.6964058, -230.49305883, -232.32509228, -230.07314787,
                -231.72602331, -228.67022594, -230.58062675, -232.05714829, -233.11207428, -231.78318343, -232.81874704,
                -231.69051811, -231.04156462, -232.13916759, -232.26584581, -232.1842805, -233.6457825, -233.16723773,
                -232.74395021, -232.81595047, -232.41401241, -231.75681963, -232.08839623, -229.31978215, -230.24534805,
                -230.86939434, -232.939754]

        self.sset6.set_property_values(property_name='total_energy_emt', property_vals=temp)
        print('energies=', temp)
        r = 3.5
        from clusterx.clusters.clusters_pool import ClustersPool
        cpool = ClustersPool(platt2, npoints=[0, 1, 2, 3, 4], radii=[0, 0, r, r, r])
        print(len(cpool), " clusters were generated.")
        from clusterx.model import ModelBuilder

        mb = ModelBuilder(selector_type="linreg", selector_opts={'clusters_sets': 'size'},
                          estimator_type="skl_LinearRegression", estimator_opts={"fit_intercept": False})
        self.cemodel1 = mb.build(self.sset6, cpool, "total_energy_emt")  # Build CE model using the training data set
        cpool_opt1 = mb.get_opt_cpool()

        # cemodel1.report_errors(sset6)
        # cpool_opt1.display_info(ecis=cemodel1.get_ecis())
        from clusterx.visualization import plot_optimization_vs_number_of_clusters
        from clusterx.visualization import plot_predictions_vs_target
        # plot_optimization_vs_number_of_clusters(mb.get_selector(),scale=0.7)
        # plot_predictions_vs_target(sset6,cemodel1,"total_energy_emt",scale=0.7)

        # Montecarlo




    def make_sample_folder(self):
        mdict={}
        for i in range(0,20):
            mc1 = MonteCarlo(self.cemodel1, self.scell2, {0: [46, 79], 1: [0, 8]}, no_of_swaps=1)
            mc2 = mc1.metropolis([8.6 * 10 ** -5, 1000], 20000, self.sset6.get_structure(0).get_atomic_numbers())
            nid = mc2.get_sampling_step_nos()
            os.mkdir('Poscar samples/{}/'.format(i))
            vasp.write_vasp('Poscar samples/{}/POSCAR'.format(i), mc2.get_structure_at_step(5000), sort=True)
            print('trajectory entries', nid)
            print('Number of accepted=', len(mc2.get_model_total_energies()))
            print(mc2.get_structure_at_step(5000).get_atoms())
            print(mc2.get_sampling_step_entry_at_step(5000)['model_total_energy'])
            mdict[i] = mc2.get_sampling_step_entry_at_step(5000)['model_total_energy']
            print(self.cemodel1.predict(mc2.get_structure_at_step(5000)))
        print(mdict)


MCsim().make_sample_folder()