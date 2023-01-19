sset3.set_calculator(EMT2())
print(sset3.calculate_property("total_energy_emt")) # Calculate energies with Effective Medium Theory calculator of ASE
g=sset3.calculate_property("total_energy_emt")
for coords in sset3:

    print(coords.get_potential_energy())
# r=3.5
# from clusterx.clusters.clusters_pool import ClustersPool
# cpool = ClustersPool(platt2, npoints=[0,1,2,3,4], radii=[0,0,r,r,r])
# print(len(cpool)," clusters were generated.")
#
#
#
#
#
#
# from clusterx.model import ModelBuilder
#
# mb = ModelBuilder(selector_type="linreg",selector_opts={'clusters_sets':'size'},estimator_type="skl_LinearRegression",estimator_opts={"fit_intercept":False})
# cemodel1 = mb.build(sset3, cpool, 'total_energy_emt') #Build CE model using the training data db
# cpool_opt1 = mb.get_opt_cpool()
#
# cemodel1.report_errors(sset3)
# cpool_opt1.display_info(ecis=cemodel1.get_ecis())
#
#
# from clusterx.visualization import plot_optimization_vs_number_of_clusters
# from clusterx.visualization import plot_predictions_vs_target
# plot_optimization_vs_number_of_clusters(mb.get_selector(),scale=0.7)
# plot_predictions_vs_target(sset3,cemodel1,"total_energy_emt",scale=0.7)
# #plot_property_vs_concentration(sset3,property_name="total_energy_emt",cemodel=cemodel1,refs=ref_en,scale=0.7)
#
#



# with open('slab.xyz', 'r') as f:
#     wordlist = [line.split(None, 1)[0] for line in f]
# wordlist.pop(0)
# wordlist.pop(0)
# print(wordlist)
#
# #iterate though list and turn element name into number pd=0 au=1
#
#
#
# wordlist = list(map(lambda x: x.replace('Pd', '0'), wordlist))
# wordlist = list(map(lambda x: x.replace('Au', '1'), wordlist))
# wordlist = list(map(int, wordlist))
#
# print(np.array(wordlist))
#
# setattr(a,'sigmas',np.array(wordlist))
# print(a.sigmas)