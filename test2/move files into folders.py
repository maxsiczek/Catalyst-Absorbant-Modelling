import shutil, os

for i in range(0,100):
    files = ['INCAR', 'KPOINTS', 'POTCAR','relax.qs']
    for f in files:
        shutil.copy(f, '{}'.format(i))