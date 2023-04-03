# # total number of iteration you want to run the code for
# N_total = 100
#
# # number of iteration after which you attempt a GCMC move
# N_GCMC = 100
#
# init_struc = ...
# GCMC_struc = init_struc
# for i in range(0,N_total):
#
#     # 1 run regular MC for N_GCMC steps
#     # from cell
#     MC(GCMC_struc,N_GCMC)
#
#     # 2 GCMC move
#     # take the last structure and
#     # randomly choose a surface site
#     # draw a integer between 1-16
#     # change the identity of the site. Meaning if its O change it to X and vice versa.
#
#     # take this new structure
#     # and assign it to GCMC_struc
#
#     GCMC_struc = struc_generated_from_add_delete_move
#
#
#     # accept/reject (later)


list=[4,8,12,16]
print(list)
list[3]=100
print(list)
list2=[4,8,12,16]
print([list,list2])