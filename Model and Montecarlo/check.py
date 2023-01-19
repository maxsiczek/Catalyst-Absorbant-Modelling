import os

nfiles=100

for i in range(0,nfiles):
    os.chdir(str(i))
    with open('OUTCAR', 'r') as f:
        Good=False
        for line in f:
            if 'reached required accuracy' in line:
                Good=True
                print(str(i)+' Good')
    if Good==False:
        #os.system('mv CONTCAR POSCAR')
        #os.system('sbatch relax.qs')
        print(str(i)+' Bad')
    os.chdir('..')

