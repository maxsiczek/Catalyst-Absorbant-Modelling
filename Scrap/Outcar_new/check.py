import os

nfiles=60
a=[]
for i in range(0,nfiles):
    os.chdir(str(i))
    with open('OUTCAR', 'r') as f:
        Good=False
        for line in f:
            if 'reached required accuracy' in line:
                Good=True
                print(str(i)+' Good')
                a.append(i)
    if Good==False:
        print(str(i)+' Bad')
    os.chdir('..')
print(a)


