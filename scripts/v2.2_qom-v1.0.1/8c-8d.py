# dependencies
import numpy as np

# qom modules
from qom.ui import init_log
from qom.ui.plotters import MPLPlotter

# initialize logger
init_log()

# frequently used variables
dim = 601
xs = np.around(np.linspace(-3, 3, dim), 2)
ys = np.around(np.linspace(-3, 3, dim), 2)
Wp = np.zeros((dim, dim), dtype=np.float_)
Wm = np.zeros((dim, dim), dtype=np.float_)
plotter_params = {
    'type'      : 'contourf',
    'x_ticks'   : [-3, -1.5, 0, 1.5, 3],
    'y_ticks'   : [-3, -1.5, 0, 1.5, 3],
    'show_cbar' : True,
    'width'     : 5.5
}

# load data
i = 0
with open('data/v2.2_qom-v1.0.1/8c-8d.txt', 'r') as datafile:
    for line in datafile:
        content = line.split(' ')
        x = float(content[0])
        y = float(content[1])
        if x in xs and y in ys:
            Wp[int(i % dim), int(i / dim)] = np.float_(content[2])
            Wm[int(i % dim), int(i / dim)] = np.float_(content[3])
            i += 1
    datafile.close()
# plus mode
plotter_params['x_label'] = '$q_{+}$'
plotter_params['y_label'] = '$p_{+}$'
plotter_params['cbar_ticks'] = [0.0, 0.1, 0.2, 0.3, 0.4]
plotter = MPLPlotter(
    axes={
        'X' : xs,
        'Y' : ys
    },
    params=plotter_params
)
plotter.update(
    vs=Wp
)
plotter.show()
# minus mode
plotter_params['x_label'] = '$q_{-}$'
plotter_params['y_label'] = '$p_{-}$'
plotter_params['cbar_ticks'] = [0.00, 0.05, 0.1, 0.15, 0.20]
plotter = MPLPlotter(
    axes={
        'X' : xs,
        'Y' : ys
    },
    params=plotter_params
)
plotter.update(
    vs=Wm
)
plotter.show()