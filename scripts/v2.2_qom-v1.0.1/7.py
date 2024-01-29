# dependencies
import numpy as np
import scipy.constants as sc
import scipy.linalg as sl

# qom modules
from qom.solvers.differential import ODESolver
from qom.ui import init_log
from qom.ui.plotters import MPLPlotter

# initialize logger
init_log()

# system parameters
omega_m = 2 * np.pi * 1e7
T       = 0.4e-3
Delta_0 = omega_m
E_0     = 6.4e3 * omega_m
g_0     = 3.5e-5 * omega_m
gamma   = 1e-5 * omega_m
kappa   = 7.5e-1 * omega_m
Omega   = 2 * omega_m
tau     = 2 * np.pi / Omega
ts      = np.linspace(0, 50 * tau, 5001)
rates   = np.zeros(4, dtype=np.float_)

# noise matrix
if T == 0:
    n = 0
else:
    n = 1 / (np.exp(sc.hbar*omega_m/sc.k/T) - 1)
D = [   [   kappa/2,   0,          0,                  0                   ],
        [   0,          kappa/2,   0,                  0                   ],
        [   0,          0,          gamma*(2*n + 1)/2, 0                   ],
        [   0,          0,          0,                  gamma*(2*n + 1)/2  ]   ]

# function for mode rates
def func_ode(t, v, c):
    # extract modes
    q, p = v[0:2]
    alpha = v[2] + 1j * v[3]
    # update amplitude
    E = E_0 + c[0] * (np.exp(-1j * Omega * t) + np.exp(1j * Omega * t))
    # update rates
    rates[0] = omega_m * p
    rates[1] = -omega_m * q - gamma / 2 * p + np.sqrt(2) * g_0 * np.real(np.conjugate(alpha) * alpha)
    dalpha_dt= -(1j * Delta_0 + kappa / 2) * alpha + 1j * np.sqrt(2) * g_0 * alpha * q + E
    rates[2] = np.real(dalpha_dt)
    rates[3] = np.imag(dalpha_dt)

    return rates

# initialize ODE solver
solver = ODESolver(
    func=func_ode,
    params={
        'show_progress' : True,
        'ode_method'    : 'vode'
    }
)
# solve ODE without modulation
E_1 = 0.0
vs_0 = solver.solve(
    T=ts,
    iv=np.zeros(4, dtype=np.float_),
    c=[E_1]
)
# solve ODE without modulation
E_1 = 3.2e3 * omega_m
vs_1 = solver.solve(
    T=ts,
    iv=np.zeros(4, dtype=np.float_),
    c=[E_1]
)

# plot optical mode
plotter = MPLPlotter(
    axes={},
    params={
        'type'      : 'lines',
        'sizes'     : [1, 2],
        'x_label'   : 'Re[$\\langle a \\rangle$]',
        'x_ticks'   : [-4000, 2000, 8000],
        'v_label'   : 'Im[$\\langle a \\rangle$]',
        'v_ticks'   : [-9000, -6000, -3000],
        'height'    : 4.0
    }
)
plotter.update(
    vs=[vs_1[:, 3], vs_1[-1000:, 3]],
    xs=[vs_1[:, 2], vs_1[-1000:, 2]]
)
plotter.show()

# plot mechanical mode
plotter = MPLPlotter(
    axes={},
    params={
        'type'      : 'lines',
        'sizes'     : [1, 2],
        'x_label'   : '$\\langle q \\rangle$',
        'x_ticks'   : [1000, 2500, 4000],
        'v_label'   : '$\\langle p \\rangle$',
        'v_ticks'   : [-2000, 0, 2000],
        'height'    : 4.0
    }
)
plotter.update(
    vs=[vs_1[:, 1], vs_1[-1000:, 1]],
    xs=[vs_1[:, 0], vs_1[-1000:, 0]]
)
plotter.show()

# function to obtain the entanglement
def get_entanglement(v):
    # update parameters
    g = g_0 * (v[2] + 1j * v[3])
    Delta = Delta_0 - np.sqrt(2) * g_0 * np.real(v[0])
    # drift matrix
    A = [   [   -kappa/2,       Delta,          -2*np.imag(g),  0           ],
            [   -Delta,         -kappa/2,       2*np.real(g),   0           ],
            [   0,              0,              -gamma/2,       omega_m     ],
            [   2*np.real(g),   2*np.imag(g),   -omega_m,       -gamma/2    ]   ]
    # invariants
    V = sl.solve_lyapunov(np.array(A), -np.array(D))
    A = [   [   V[0][0],    V[0][1]     ],
            [   V[1][0],    V[1][1]     ]   ]
    B = [   [   V[2][2],    V[2][3]     ],
            [   V[3][2],    V[3][3]     ]   ]
    C = [   [   V[0][2],    V[0][3]     ],
            [   V[1][2],    V[1][3]     ]   ]
    Sigma   = np.linalg.det(A) + np.linalg.det(B) - 2*np.linalg.det(C)
    # return entanglement
    return np.max([0, - np.log(2*np.sqrt((Sigma - np.sqrt(Sigma**2 - 4*np.linalg.det(V)))/2))])

# calculate entanglement
E_N_0 = [get_entanglement(vs_0[i]) for i in range(4500, 5001)]
E_N_1 = [get_entanglement(vs_1[i]) for i in range(4500, 5001)]

# plot entanglement
plotter = MPLPlotter(
    axes={},
    params={
        'type'          : 'lines',
        'x_label'       : '$t / \\tau$',
        'x_tick_labels' : [45, 46, 47, 48, 49, 50],
        'x_ticks'       : ts[4500:5001:100],
        'v_label'       : '$E_{N}$',
        'v_ticks'       : [0.1, 0.2, 0.3],
        'width'         : 10.0,
        'height'        : 4.0
    }
)
plotter.update(
    vs=[E_N_0, E_N_1],
    xs=ts[4500:5001]
)
plotter.show()