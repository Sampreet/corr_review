from matplotlib import colors as mc
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import scipy.linalg as sl

# plt.figure()
# plt.rcParams.update({'font.size': 18, 'axes.linewidth': 2})
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
# plt.xlabel(r'$\Delta / \omega_m$', fontsize=18)
# plt.ylabel(r'$E_N$', fontsize=18)

omega_m     = 2 * np.pi * 1e7     
kappa       = 7.5e-1 * omega_m
gamma       = 1e-5 * omega_m
T           = 0.4e-3
# g_0s        = np.array([1.5e3, 2.5e3])
g_0s        = np.array([3e-5, 4e-5]) * omega_m
E_0         = 6.4e3 * omega_m
Deltas      = np.linspace(0, 3, 1001) * omega_m
colors      = ['g', 'b']
linestyles  = ['--', '-.']

E_Ns = []

for i in range(len(g_0s)):
    g_0 = g_0s[i]
    E_N = []
    for Delta in np.nditer(Deltas):
        alpha = np.abs(E_0 / (kappa / 2 + 1j * Delta))
        g = g_0 * alpha

        expr_1 = Delta**4*gamma*kappa + Delta**2*gamma**3*kappa/2 + Delta**2*gamma**2*kappa**2 + Delta**2*gamma*kappa**3/2 - 2*Delta**2*gamma*kappa*omega_m**2 + 4*Delta*g**2*gamma**2*omega_m + 8*Delta*g**2*gamma*kappa*omega_m + 4*Delta*g**2*kappa**2*omega_m + gamma**5*kappa/16 + gamma**4*kappa**2/4 + 3*gamma**3*kappa**3/8 + gamma**3*kappa*omega_m**2/2 + gamma**2*kappa**4/4 + gamma**2*kappa**2*omega_m**2 + gamma*kappa**5/16 + gamma*kappa**3*omega_m**2/2 + gamma*kappa*omega_m**4

        expr_2 = 4*Delta**2*gamma**2 + 16*Delta**2*omega_m**2 - 64*Delta*g**2*omega_m + gamma**2*kappa**2 + 4*kappa**2*omega_m**2

        if (expr_1 > 0 and expr_2 > 0):
            A = [   [   -kappa/2,   Delta,      0,          0           ],
                    [   -Delta,     -kappa/2,   2*g,        0           ],
                    [   0,          0,          -gamma/2,   omega_m     ],
                    [   2*g,        0,          -omega_m,   -gamma/2    ]   ]
            if T == 0:
                n = 0
            else:
                n = 1 / (np.exp(sc.hbar*omega_m/sc.k/T) - 1)
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

            if Delta / omega_m == 0.855:
                print(en)
        else:
            E_N.append(0)

    E_Ns.append(E_N)

    # plt.plot(Deltas / omega_m, E_N, color=colors[i], linestyle=linestyles[i], linewidth=2, label=r'$g / \omega_m = {:1.0f}'.format(g_0s[i] / omega_m * 1e5) + '\\times 10^{-5}$')
# print((Deltas/omega_m).tolist())
# plt.legend(frameon=False)
# plt.show()

from qom.ui import figure
from qom.utils import axis
import seaborn as sns

plot_params = {
    'y_label': '$E_{N}$',
    'type': 'lines',
    'label_font_size': 20,
    'tick_font_size': 16
}

colors = sns.diverging_palette(250, 15, s=75, l=40, n=11, center='light', as_cmap=False)

x_dict = {
    'var': 'xs',
    'label': '$\\Delta / \\omega_{m}$',
    'values': np.around(Deltas / omega_m, 3),
    'tick_labels': [0, 1, 2, 3]
}
X = axis.StaticAxis(x_dict)

z_dict = {
    'var': 'zs',
    'label': '$g / \\omega_m$',
    'unit': '$\\times 10^{-5}$',
    'values': [3, 4],
    'colors': [colors[0], colors[-1]],
    'linestyles': ['-', '-']
}
Z = axis.StaticAxis(z_dict)

plotter = figure.Plotter(plot_params, X, Z=Z)

X_v = axis.DynamicAxis([2, len(Deltas)])
X_v.values[0] = X.values
X_v.values[1] = X.values

V = axis.DynamicAxis([2, len(Deltas)])
V.values = E_Ns

plotter.update(X=X_v, Y=V)