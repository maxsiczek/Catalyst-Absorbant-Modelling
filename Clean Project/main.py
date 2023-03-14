import os

from Catalyst_Analysis import TsetAdsorbate
c=TsetAdsorbate(['Pd','Au'],100,None,'O',8,'fcc')
# os.mkdir('test')
os.chdir('test')
c.sset_obj()
# c.create_poscar_files('POTCAR','structure_oxygen.pkl')
