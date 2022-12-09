import sys
import numpy as np
import matplotlib.pyplot as plt

names = sys.argv[1:]
data = np.loadtxt(names[0], delimiter='\t', unpack=True)
x = np.arange(1, len(data[0]) + 1, 1, dtype=int)

plt.style.use('science')
fig, ax = plt.subplots()
plt.title(names[2])
plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0.0))

ax.errorbar(x, data[0], yerr=data[1], fmt='.', color='blue')

ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
if len(x) > 10:
    ax.set_xticks(np.arange(1, len(x) + 1, len(x) // 10))
else:
    ax.set_xticks(np.arange(1, len(x) + 1, 1))

ax.set_xlabel(names[3])
ax.set_ylabel(names[4])
ax.legend(loc='upper left')

plt.tight_layout()
plt.savefig(names[1], dpi=250)

# Path: tools/plots/measurements_error_bars.py
# Usage: python3 measurements_error_bars.py data_path/.txt image_path/output.png 'title' 'x-axis' 'y-axis'