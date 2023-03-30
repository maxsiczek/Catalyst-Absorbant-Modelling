#Extracts the energies from OUTCAR files in the same folder and makes them into a list. Number of OUTCAR files must be specified.


from ase.io import vasp
import os

nfiles=100
Goodlist=[]
for i in range(0,nfiles):
    os.chdir(str(i))
    with open('OUTCAR', 'r') as f:
        Good=False
        for line in f:
            if 'reached required accuracy' in line:
                Good=True
                print(str(i)+' Good')
                Goodlist.append(i)
    if Good==False:
        print(str(i)+' Bad')
    os.chdir('..')
print(Goodlist)



Elist= []
for i in Goodlist:
    os.chdir(str({}).format(i))
    coords=vasp.read_vasp_out('OUTCAR')
    os.chdir('..')
    # print(atoms_object.get_cell())
    Elist.append(coords.get_potential_energy())
print(Elist)
print(str(len(Elist))+" energies")
