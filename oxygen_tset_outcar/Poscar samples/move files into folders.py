#pushes all required files into the folders
import shutil, os

for i in range(0,20):
    files = ['INCAR', 'KPOINTS', 'POTCAR','relax.qs']
    for f in files:
        shutil.copy(f, '{}'.format(i))