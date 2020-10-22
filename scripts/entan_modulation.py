from matplotlib import colors as mc, cm
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import scipy.integrate as si
import scipy.linalg as sl
import seaborn as sns

def model (t, v, c):
    omega_m, Delta_0, gamma, kappa, g_0, E_0, E_1, Omega = c
    rates = []
    E = E_0 + E_1 * (np.exp(-1j * Omega * t) + np.exp(1j * Omega * t))
    rates.append(omega_m * v[1])
    rates.append(-omega_m * v[0] - gamma / 2 * v[1] + np.sqrt(2) * g_0 * v[2] * np.conjugate(v[2]))
    rates.append(-(1j * Delta_0 + kappa / 2) * v[2] + 1j * np.sqrt(2) * g_0 * v[2] * v[0] + E)

    return rates

omega_m = 2 * np.pi * 1e7
Delta_0 = omega_m
kappa   = 7.5e-1 * omega_m
gamma   = 1e-5 * omega_m
T       = 0.4e-3
# g_0     = 2e3
g_0     = 3.5e-5 * omega_m
E_0     = 6.4e3 * omega_m
E_1     = 0.5 * E_0
if T == 0:
    n = 0
else:
    n = 1 / (np.exp(sc.hbar*omega_m/sc.k/T) - 1)
Omega   = 2 * omega_m
tau     = 2 * np.pi / Omega
teval   = np.linspace(0, 50 * tau, 5001)

t = 0
v = [0, 0, 0]
c = [omega_m, Delta_0, gamma, kappa, g_0, E_0, E_1, Omega]
integrator = si.ode(model)
integrator.set_integrator('zvode')
integrator.set_initial_value(v, t)
integrator.set_f_params(c)

V = [   [1/2,   0,      0,      0       ],
        [0 ,    1/2,    0,      0       ],
        [0 ,    0,      n + 1/2,0       ],
        [0 ,    0,      0,      n + 1/2 ]   ]
V2I = []
V0R = []
V1R = []
E_N = [0]
for i in range(1, len(teval)):
    t = teval[i]
    v = integrator.integrate(t)
    a = v[2]
    g = g_0 * a
    Delta = Delta_0 - np.sqrt(2) * g_0 * np.real(v[0])
    A = [   [   -kappa/2,       Delta,          -2*np.imag(g),  0           ],
            [   -Delta,         -kappa/2,       2*np.real(g),   0           ],
            [   0,              0,              -gamma/2,       omega_m     ],
            [   2*np.real(g),   2*np.imag(g),   -omega_m,       -gamma/2    ]   ]
    D = [   [   kappa/2,   0,          0,                  0                   ],
            [   0,          kappa/2,   0,                  0                   ],
            [   0,          0,          gamma*(2*n + 1)/2, 0                   ],
            [   0,          0,          0,                  gamma*(2*n + 1)/2  ]   ]
    V = sl.solve_lyapunov(np.array(A), -np.array(D))
    A = [   [   V[0][0],    V[0][1]     ],
            [   V[1][0],    V[1][1]     ]   ]
    B = [   [   V[2][2],    V[2][3]     ],
            [   V[3][2],    V[3][3]     ]   ]
    C = [   [   V[0][2],    V[0][3]     ],
            [   V[1][2],    V[1][3]     ]   ]
    Sigma   = np.linalg.det(A) + np.linalg.det(B) - 2*np.linalg.det(C)
    eta_m   = np.sqrt((Sigma - np.sqrt(Sigma**2 - 4*np.linalg.det(V)))/2)
    en      = np.max([0, - np.log(2*eta_m)])
    E_N.append(en)

dim = 501
E_N0 = np.ones(dim) * 0.18966077544510937

# plt.figure()
# plt.rcParams.update({'font.size': 16, 'axes.linewidth': 2})
# plt.plot(teval[-1000:] / tau, E_N0, color='k', linestyle='--', linewidth=2)
# plt.plot(teval[-1000:] / tau, E_N[-1000:], color='r', linestyle='-', linewidth=2)
# plt.xlabel(r'$t / \tau$', fontsize=24)
# plt.ylabel(r'$E_N$', fontsize=24)
# # plt.legend(['witout modulation', 'with modulation'], loc='upper right')

# plt.show()

from qom.ui import figure
from qom.utils import axis
plot_params = {
    'y_label': '$E_{N}$',
    'y_ticks': [0.1, 0.2, 0.3],
    'type': 'lines',
    'show_legend': False,
    'label_font_size': 20,
    'tick_font_size': 16
}

colors = sns.diverging_palette(250, 15, s=75, l=40, n=11, center='light', as_cmap=False)
print(colors[-1])

x_dict = {
    'var': 'xs',
    'label': '$t / \\tau$',
    'values': np.around(teval[- dim :] / tau, 3)
}
X = axis.StaticAxis(x_dict)

z_dict = {
    'var': 'zs',
    'label': 'modulation',
    'values': [0, 1],
    'colors': [colors[0], colors[-1]],
    'linestyles': ['--', '-']
}
Z = axis.StaticAxis(z_dict)

plotter = figure.Plotter(plot_params, X, Z=Z)

X_v = axis.DynamicAxis([2, dim])
X_v.values[0] = X.values
X_v.values[1] = X.values

V = axis.DynamicAxis([2, dim])
V.values[0] = E_N0
V.values[1] = E_N[- dim :]

plotter.update(X=X_v, Y=V)