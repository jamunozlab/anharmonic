# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from ase import Atoms
from ase.visualize import view
from ase.io import write
from copy import copy
from scipy.spatial.distance import cdist
import pandas as pd

d = 14.07
    
seed = [float(x)/8 for x in range(0, 8, 2)]
coords = []
for i in seed:
    for j in seed:
        for k in seed:
            coords.append((i, j, k))
            coords.append((i+.125, j+.125, k+.125))

coords = np.array(coords)
coords = coords*d

atoms = Atoms(['Zr' for x in range(len(coords))], positions=coords, cell=[d, d, d], pbc=True)

#view(atoms)
#write('-', atoms, format='espresso-in')
#write('/Users/jamunoz/OneDrive - University of Texas at El Paso/git/anharmonic/bcc_Zr_128.cif', atoms, format='cif')

init_pos = copy(atoms.positions)

#print(len(cdist(atoms.positions, init_pos)))

#write('-', atoms, format='cif')

dists = []
filenames = []
path = '/Users/jamunoz/OneDrive - University of Texas at El Paso/git/anharmonic/bcc_Zr_128_rattled/'
for i in range(2):
    atoms = Atoms(['Zr' for x in range(len(coords))], positions=coords, cell=[d, d, d], pbc=True)
    atoms.rattle(seed=i)
    filename = 'bcc_Zr_128_rattled_' + str(i) + '.in'
    write(path + filename, atoms, format='espresso-in')
    dists.append( np.sum(cdist(atoms.positions, init_pos)) / 1.512e5)
    filenames.append(filename)

energy = pd.Series(dists) - 1
energy =  energy * 100

aDict = {'filename': filenames, 'energy': energy}
df = pd.DataFrame(aDict)

path = '/Users/jamunoz/OneDrive - University of Texas at El Paso/git/anharmonic/'
df.to_csv(path + 'bcc_Zr_128_rattled_energies.csv')

print(df)
