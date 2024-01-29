# dependencies
import numpy as np
import scipy.constants as sc
import scipy.linalg as sl

# qom modules
from qom.ui.plotters import MPLPlotter

# system parameters
omega_m = 2*np.pi*1e7  
T       = 0.4e-3
Delta   = omega_m
E_0     = 6.4e3 * omega_m
gamma   = 1e-5 * omega_m              
kappa   = 7.5e-1 * omega_m
gs      = np.array([3, 4]) * 1e-5 * omega_m

# noise matrix
if T == 0:
    n = 0
else:
    n = 1 / (np.exp(sc.hbar*omega_m/sc.k/T) - 1)
D = [   [   kappa/2,   0,          0,                  0                   ],
        [   0,          kappa/2,   0,                  0                   ],
        [   0,          0,          gamma*(2*n + 1)/2, 0                   ],
        [   0,          0,          0,                  gamma*(2*n + 1)/2  ]   ]

# obtain elements
Vs = np.zeros((2, 10), dtype=np.float_)
idxs = np.triu_indices(4)
for i in range(2):
    # update parameters
    G = gs[i] * E_0 / np.sqrt(kappa**2 / 4 + Delta**2)
    # drift matrix
    A = [   [   -kappa/2,   Delta,      0,          0           ],
            [   -Delta,     -kappa/2,   2*G,        0           ],
            [   0,          0,          -gamma/2,   omega_m     ],
            [   2*G,        0,          -omega_m,   -gamma/2    ]   ]
    # covariance matrix
    _V = sl.solve_lyapunov(np.array(A), -1 * np.array(D))
    Vs[i] = _V[idxs]
    Vs[i][6] *= -1 

# plotter
plotter = MPLPlotter(
    axes={},
    params={
        'type'          : 'scatters',
        'sizes'         : [20, 20],
        'styles'        : ['x', 'o'],
        'x_label'       : '',
        'x_limits'      : [-0.5, 9.5],
        'x_tick_labels' : ['$V_{' + str(idx) + '}$' if idx != 24 else '$-V_{24}$' for idx in ((idxs[0] + 1) * 10 + idxs[1] + 1)],
        'x_ticks'       : list(range(10)),
        'v_label'       : '',
        'v_limits'      : [-0.1, 0.7],
        'v_ticks'       : [0.0, 0.2, 0.4, 0.6],
        'width'         : 5.0,
        'height'        : 2.0
    }
)
plotter.update(
    vs=Vs,
    xs=list(range(10))
)
plotter.show()
