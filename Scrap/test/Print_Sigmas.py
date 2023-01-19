#Prints the sigmas and element list from each folder containing a xyz file


import os
import numpy as np

nfiles=5
for i in range(0,nfiles):

    os.chdir(str(i))
    with open('{}.xyz'.format(i), 'r') as f:
        wordlist = [line.split(None, 1)[0] for line in f]
    wordlist.pop(0)
    wordlist.pop(0)
    print(wordlist)


    #iterate though list and turn element name into number pd=0 au=1
    wordlist = list(map(lambda x: x.replace('Pd', '0'), wordlist))
    wordlist = list(map(lambda x: x.replace('Au', '1'), wordlist))
    wordlist = list(map(int, wordlist))

    print(np.array(wordlist))
    os.chdir('..')