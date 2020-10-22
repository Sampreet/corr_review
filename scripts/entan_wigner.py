from matplotlib import colors as mc, cm
import matplotlib.pyplot as plt
import numpy as np

dim = 601
x_data = np.around(np.linspace(-3, 3, dim), 2)
y_data = np.around(np.linspace(-3, 3, dim), 2)
Wp = []
Wm = []

with open('data/wig_0.01_300.45.txt', 'r') as datafile:
    for line in datafile:
        content = line.split(' ')
        x = float(content[0])
        y = float(content[1])
        if x in x_data and y in y_data:
            Wp.append(float(content[2]))
            Wm.append(float(content[3]))
    datafile.close()
X, Y = np.meshgrid(x_data, y_data)
Wp = np.transpose(np.array(Wp).reshape(dim, dim))
Wm = np.transpose(np.array(Wm).reshape(dim, dim))
colors = ['b', 'w', 'r']
cmap = mc.LinearSegmentedColormap.from_list('hello', colors, 101)

# fig = plt.figure()
# axes = fig.gca()
# plt.rcParams.update({'font.size': 16})
# plt.xticks(np.linspace(-5, 5, 5), fontsize=16)
# plt.yticks(np.linspace(-5, 5, 5),fontsize=16)
# axes.set_xlabel(r'$\delta q_{+}$', fontsize=14)
# axes.set_ylabel(r'$\delta p_{+}$', fontsize=14)
# fig_p = axes.contourf(X, Y, Wp, cmap=cmap)
# fig.colorbar(fig_p, shrink=1, aspect=15)

# fig = plt.figure()
# axes = fig.gca()
# plt.rcParams.update({'font.size': 16})
# plt.xticks(np.linspace(-5, 5, 5), fontsize=16)
# plt.yticks(np.linspace(-5, 5, 5),fontsize=16)
# axes.set_xlabel(r'$\delta q_{-}$', fontsize=14)
# axes.set_ylabel(r'$\delta p_{-}$', fontsize=14)
# fig_m = axes.contourf(X, Y, Wm, cmap=cmap)
# fig.colorbar(fig_m, shrink=1, aspect=15)

# plt.show()

from qom.ui import figure
from qom.utils import axis
plot_params = {
    'cbar_label': '$W (q_{+}, p_{+})$',
    'cbar_ticks': [0, 1, 2, 3, 4],
    'type': 'contourf',
    'label_font_size': 24,
    'tick_font_size': 20,
    'show_cbar': True
}

x_dict = {
    'var': 'xs',
    'label': '$q_{+}$',
    'values': x_data
}
X = axis.StaticAxis(x_dict)

y_dict = {
    'var': 'ys',
    'label': '$p_{+}$',
    'values': y_data
}
Y = axis.StaticAxis(y_dict)

plotter = figure.Plotter(plot_params, X, Y)

Z = axis.DynamicAxis([len(y_data), len(x_data)])
Z.values = Wp

plotter.update(Z=Z)