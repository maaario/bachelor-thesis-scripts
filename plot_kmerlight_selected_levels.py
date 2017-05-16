# Plot levels selected by kmerlight

import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({'font.size': 12})

levels = []
abundances = []

with open('kmerlight_levels.txt', 'r') as f:
    for line in f.readlines():
        i, w, ti, fi = [int(float(x)) for x in line.split()]

        abundances.append(i)
        levels.append(w)

plt.scatter(abundances, levels, alpha=0.2)
plt.xlabel('i')
plt.ylabel('w')
plt.ylim([6.5, 12.5])
plt.show()
