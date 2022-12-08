import sys
import numpy as np
import matplotlib.pyplot as plt

names = sys.argv[1:]
print(names)
x_data = np.loadtxt(names[0], delimiter='\t', unpack=True)
y_data = np.loadtxt(names[1], delimiter='\t', unpack=True)

plt.style.use('science')
fig, ax = plt.subplots()
plt.title(names[3])
plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0.0))

ax.errorbar(x_data[0], y_data[0], xerr=x_data[1], yerr=y_data[1], fmt='.', color='blue')

ax.set_xlabel(names[4])
ax.set_ylabel(names[5])
ax.legend(loc='upper left')

plt.tight_layout()
plt.savefig(names[2], dpi=250)
