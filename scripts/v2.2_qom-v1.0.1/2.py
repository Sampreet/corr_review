# dependencies
import numpy as np

# qom modules
from qom.solvers.measure import get_Wigner_distributions_single_mode
from qom.ui import init_log
from qom.ui.plotters import MPLPlotter

# initialize logger
init_log()

# frequently used variables
qs = np.linspace(-5, 5, 401)
ps = np.linspace(-5, 5, 401)
plotter_params = {
    'type'      : 'surface_cz',
    'x_label'   : '$q$',
    'x_limits'  : [-5.5, 5.5],
    'x_tick_pad': 0,
    'x_ticks'   : [-5, -2.5, 0, 2.5, 5],
    'y_label'   : '$p$',
    'y_limits'  : [-5.5, 5.5],
    'y_tick_pad': 0,
    'y_ticks'   : [-5, -2.5, 0, 2.5, 5],
    'v_label'   : '$W(q, p)$',
    'v_ticks'   : [0.0, 0.1, 0.2, 0.3]
}

# vacuum state
plotter = MPLPlotter(
    axes={
        'X' : qs,
        'Y' : ps
    },
    params=plotter_params
)
plotter.update(
    vs=get_Wigner_distributions_single_mode(
        Corrs = np.array([[[1/2, 0], [0, 1/2]]]),
        params={
            'show_progress' : True,
            'wigner_xs'     : qs,
            'wigner_ys'     : ps
        }
    )[0, 0]
)
plotter.show()

# coherent state
alpha = 1 + 1j
plotter = MPLPlotter(
    axes={
        'X' : qs,
        'Y' : ps
    },
    params=plotter_params
)
plotter.update(
    vs=get_Wigner_distributions_single_mode(
        Corrs = np.array([[[1/2, 0], [0, 1/2]]]),
        params={
            'show_progress' : True,
            'wigner_xs'     : qs - np.real(alpha),
            'wigner_ys'     : ps - np.imag(alpha)
        }
    )[0, 0]
)
plotter.show()

# squeezed state
r = 0.5
V_s = [[np.exp(-2*r) / 2, 0], [0, np.exp(2*r) / 2]]
# plotter
plotter = MPLPlotter(
    axes={
        'X' : qs,
        'Y' : ps
    },
    params=plotter_params
)
plotter.update(
    vs=get_Wigner_distributions_single_mode(
        Corrs = np.array([V_s]),
        params={
            'show_progress' : True,
            'wigner_xs'     : qs,
            'wigner_ys'     : ps
        }
    )[0, 0]
)
plotter.show()