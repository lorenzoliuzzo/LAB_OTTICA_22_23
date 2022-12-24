import sys
import numpy as np
import matplotlib.pyplot as plt

names = sys.argv[1:]
data = np.loadtxt(names[0], delimiter='\t', unpack=True)
x = np.linspace(1, len(data[0]), len(data[0]), dtype=int)

plt.style.use('science')
fig, ax = plt.subplots()

plt.title(names[2])
plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0.0))

media = np.mean(data[0])
std = np.std(data[0])
ax.errorbar(x, data[0], yerr=data[1], fmt='.', color='blue')
ax.axhline(media, color='red', label='$n_m$')
ax.axhline(media + std, linestyle='dashed', linewidth=1, color='red', label='$n_m \pm \sigma$')
ax.axhline(media - std, linestyle='dashed', linewidth=1, color='red')

ax.set_xticks(x)
ax.set_xlabel(names[3])
ax.set_ylabel(names[4])
ax.legend(loc='best')

plt.tight_layout()
plt.savefig(names[1], dpi=250)

# Usage: python3 ../../tools/plots/compare_results_vertical.py data/.dat images/.png 'title' 'x-axis' 'y-axis'