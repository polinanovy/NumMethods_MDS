import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.cm
from matplotlib import colors
import matplotlib
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams.update({'font.size': 18})
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Times New Roman'
matplotlib.rcParams['mathtext.it'] = 'Times New Roman:italic'
matplotlib.rcParams['mathtext.bf'] = 'Times New Roman:bold'

with open('RESULT') as table:
	line1 = table.readline()
	N = int(line1.split()[0])
	M= int(line1.split()[1])
	dt = float(line1.split()[2])
	t_stop = float(line1.split()[3])
	u = np.empty([M, N])
	for index, line in enumerate(table):
		u[index,:] = line.split()

mmax = u.max()
mmin = u.min()
plt.figure(figsize=(8,6))
plt.title(r'$T~(t_{stop}=$'+f'${t_stop}~$s, '+r'$~x,~y)$')
plt.xlabel(r'$x,$ m')
plt.ylabel(r'$y,$ m')
#plt.gca().xaxis.set_major_locator(MultipleLocator(10))
#plt.gca().xaxis.set_major_formatter('{x:.0f}')
cax = plt.imshow(u, interpolation = 'none', origin = 'lower', aspect = 'equal', extent = (0,0.5,0,0.5), cmap = 'Reds')
#myticks = [int(round(mmin, -1) + 10 * i) for i in range(int((mmax - mmin) // 10 + 1))]
plt.colorbar(cax)#, ticks=myticks).ax.set_yticklabels(myticks)
plt.tight_layout()
plt.savefig('map.png', format = 'png', dpi=600)


plt.show()
