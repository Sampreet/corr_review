# dependencies
import numpy as np
import scipy.constants as sc
import scipy.linalg as sl

# qom modules
from qom.utils.loopers import run_loopers_in_parallel, wrap_looper

# system parameters
omega_m = 2*np.pi*1e7
T       = 0.4e-3
E_0     = 6.4e3 * omega_m
gamma   = 1e-5 * omega_m              
kappa   = 7.5e-1 * omega_m
gs      = np.linspace(2, 5, 901) * 1e-5 * omega_m
Deltas  = np.linspace(0, 2, 1001) * omega_m

# noise matrix
if T == 0:
    n = 0
else:
    n = 1 / (np.exp(sc.hbar*omega_m/sc.k/T) - 1)
D = [   [   kappa/2,   0,          0,                  0                   ],
        [   0,          kappa/2,   0,                  0                   ],
        [   0,          0,          gamma*(2*n + 1)/2, 0                   ],
        [   0,          0,          0,                  gamma*(2*n + 1)/2  ]   ]

# function to obtain the entanglement if the system is stable
def get_entanglement_if_stable(system_params):
    # update parameters
    Delta = system_params['Delta']
    G = system_params['g'] * E_0 / np.sqrt(kappa**2 / 4 + Delta**2)

    # calculate stability condition
    expr_1 = Delta**4*gamma*kappa + Delta**2*gamma**3*kappa/2 + Delta**2*gamma**2*kappa**2 + Delta**2*gamma*kappa**3/2 - 2*Delta**2*gamma*kappa*omega_m**2 + 4*Delta*G**2*gamma**2*omega_m + 8*Delta*G**2*gamma*kappa*omega_m + 4*Delta*G**2*kappa**2*omega_m + gamma**5*kappa/16 + gamma**4*kappa**2/4 + 3*gamma**3*kappa**3/8 + gamma**3*kappa*omega_m**2/2 + gamma**2*kappa**4/4 + gamma**2*kappa**2*omega_m**2 + gamma*kappa**5/16 + gamma*kappa**3*omega_m**2/2 + gamma*kappa*omega_m**4
    expr_2 = 4*Delta**2*gamma**2 + 16*Delta**2*omega_m**2 - 64*Delta*G**2*omega_m + gamma**2*kappa**2 + 4*kappa**2*omega_m**2

    # calculate entanglement if stable
    if (expr_1 > 0 and expr_2 > 0):
        # drift matrix
        A = [   [   -kappa/2,   Delta,      0,          0           ],
                [   -Delta,     -kappa/2,   2*G,        0           ],
                [   0,          0,          -gamma/2,   omega_m     ],
                [   2*G,        0,          -omega_m,   -gamma/2    ]   ]
        # invariants
        V = sl.solve_lyapunov(np.array(A), -1 * np.array(D))
        A = [   [   V[0][0],    V[0][1]     ],
                [   V[1][0],    V[1][1]     ]   ]
        B = [   [   V[2][2],    V[2][3]     ],
                [   V[3][2],    V[3][3]     ]   ]
        C = [   [   V[0][2],    V[0][3]     ],
                [   V[1][2],    V[1][3]     ]   ]
        Sigma   = np.linalg.det(A) + np.linalg.det(B) - 2*np.linalg.det(C)
        # entanglement
        return np.max([0, - np.log(2*np.sqrt((Sigma - np.sqrt(Sigma**2 - 4*np.linalg.det(V)))/2))])
    else:
        return np.NaN
    
if __name__ == '__main__':
    # loop and plot
    looper = run_loopers_in_parallel(
        looper_name='XYLooper',
        func=get_entanglement_if_stable,
        params={
            'show_progress'     : True,
            'X'                 : {
                'var'   : 'Delta',
                'val'   : Deltas
            },
            'Y'                 : {
                'var'   : 'g',
                'val'   : gs
            }
        },
        params_system={},
        plot=True,
        params_plotter={
            'type'              : 'contourf',
            'x_label'           : '$\\Delta / \\omega_{m}$',
            'x_tick_labels'     : [0.0, 0.5, 1.0, 1.5, 2.0],
            'x_tick_position'   : 'both-out',
            'x_ticks'           : Deltas[::250],
            'y_label'           : '$10^{5} \\times g / \\omega_{m}$',
            'y_tick_labels'     : [2, 3, 4, 5],
            'y_ticks'           : gs[::300],
            'y_tick_position'   : 'both-out',
            'show_cbar'         : True,
            'cbar_ticks'        : [0.0, 0.1, 0.2, 0.3, 0.4, 0.5],
            'width'             : 5.5,
            'annotations'       : [{
                'text'          : 'unstable',
                'xy'            : [0.21, 0.68],
                'orientation'   : 'vertical',
                'color'         : 'w'
            }]
        }
    )