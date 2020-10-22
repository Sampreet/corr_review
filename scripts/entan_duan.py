from matplotlib import colors as mc
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import scipy.linalg as sl

plt.figure()
plt.rcParams.update({'font.size': 18})
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel(r'$\Delta / \omega_m$', fontsize=18)
plt.ylabel(r'$\Delta_{EPR}$', fontsize=18)

omega_m     = 2 * np.pi * 1e7     
kappa       = 7.5e-1 * omega_m
gamma       = 1e-5 * omega_m
T           = 0.4e-3
# g_0s        = np.array([1.5e3, 2.5e3])
g_0s        = np.array([3e-5, 4e-5]) * omega_m
E_0         = 4e11
Deltas      = np.linspace(0, 3, 1001) * omega_m
colors      = ['g', 'b']
linestyles  = ['--', '-.']

g_0 = 3.5e-5 * omega_m
Duan = []
for Delta in np.nditer(Deltas):
    alpha = np.abs(E_0 / (kappa / 2 + 1j * Delta))
    g = g_0 * alpha

    expr_1 = Delta**4*gamma*kappa + Delta**2*gamma**3*kappa/2 + Delta**2*gamma**2*kappa**2 + Delta**2*gamma*kappa**3/2 - 2*Delta**2*gamma*kappa*omega_m**2 + 4*Delta*g**2*gamma**2*omega_m + 8*Delta*g**2*gamma*kappa*omega_m + 4*Delta*g**2*kappa**2*omega_m + gamma**5*kappa/16 + gamma**4*kappa**2/4 + 3*gamma**3*kappa**3/8 + gamma**3*kappa*omega_m**2/2 + gamma**2*kappa**4/4 + gamma**2*kappa**2*omega_m**2 + gamma*kappa**5/16 + gamma*kappa**3*omega_m**2/2 + gamma*kappa*omega_m**4

    expr_2 = 4*Delta**2*gamma**2 + 16*Delta**2*omega_m**2 - 64*Delta*g**2*omega_m + gamma**2*kappa**2 + 4*kappa**2*omega_m**2

    if (expr_1 > 0 and expr_2 > 0 and not Delta == 0.0):
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
        duan    = V[0][0] + V[1][1] + V[2][2] + V[3][3] + 2 * (V[0][1] - V[2][3])
        print(V)
        Duan.append(duan)
    else:
        Duan.append(0)

plt.plot(Deltas[-800:] / omega_m, Duan[-800:], color='g', linestyle='-.')
plt.show()