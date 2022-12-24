import sys
import numpy as np
import matplotlib.pyplot as plt

names = sys.argv[1:]
data = np.loadtxt(names[0], delimiter='\t', unpack=True)
y = np.linspace(1, len(data[0]), len(data[0]), dtype=int)
media = np.mean(data[0])
std = np.std(data[0])

plt.style.use('science')
plt.title(names[2])
plt.ticklabel_format(axis="x", style="sci")

fig, ax = plt.subplots()
ax.errorbar(data[0], y, xerr=data[1], fmt='.', color='blue')
ax.axvline(media, color='red', label='$\lambda_m$')
ax.axvline(media + std, linestyle='dashed', linewidth=1, color='red', label='$\lambda_m \pm \sigma$')
ax.axvline(media - std, linestyle='dashed', linewidth=1, color='red')

ax.set_yticks(y)
ax.set_xlabel(names[3])
ax.set_ylabel(names[4])
ax.legend()

plt.tight_layout()
plt.savefig(names[1], dpi=250)

# Usage: python3 ../../tools/plots/compare_results_vertical.py data/.dat images/.png 'title' 'x-axis' 'y-axis' 