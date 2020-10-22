from matplotlib import colors as mc, cm
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import scipy.linalg as sl
import seaborn as sns

omega_m     = 2*np.pi*1e7                   
kappa       = 7.5e-1 * omega_m
gamma       = 1e-5 * omega_m
Delta       = omega_m
T           = 0.4e-3
E_0         = 6.4e3 * omega_m
g_0s        = np.linspace(2, 5, 1001) * 1e-5 * omega_m
Deltas      = np.linspace(0, 2, 1001) * omega_m

_maxi = []
_mini = []
E_N = []
Ins = []
for i in range(len(g_0s)):
    print(i/len(g_0s))
    g_0 = g_0s[i]
    temp_en = []
    temp_in = []
    for Delta in np.nditer(Deltas):
        alpha = E_0 / np.sqrt(kappa**2 / 4 + Delta**2)
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
            V = sl.solve_lyapunov(np.array(A), -1 * np.array(D))
            A = [   [   V[0][0],    V[0][1]     ],
                    [   V[1][0],    V[1][1]     ]   ]
            B = [   [   V[2][2],    V[2][3]     ],
                    [   V[3][2],    V[3][3]     ]   ]
            C = [   [   V[0][2],    V[0][3]     ],
                    [   V[1][2],    V[1][3]     ]   ]
            Sigma   = np.linalg.det(A) + np.linalg.det(B) - 2*np.linalg.det(C)
            eta_m   = np.sqrt((Sigma - np.sqrt(Sigma**2 - 4*np.linalg.det(V)))/2)
            en      = np.max([0, - np.log(2*eta_m)])
            temp_en.append(en)
            temp_in.append(np.NaN)
        else:
            temp_en.append(0)
            temp_in.append(0)

    _mini.append(min(temp_en))
    _maxi.append(max(temp_en))
    E_N.append(temp_en)
    Ins.append(temp_in)
    
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams.update({
    'font.family': 'Times New Roman',
    'font.style': 'normal',
    'font.variant': 'normal',
    'font.weight': 500,
    'font.stretch': 500,
    'font.size': 20
})

fig = plt.figure()
# axes = fig.gca(projection='3d')
axes = fig.gca()
plt.rcParams.update({'font.size': 16})
plt.xticks(np.linspace(0, 2, 5), fontsize=16)
plt.yticks(np.linspace(2, 5, 4), fontsize=16)
axes.set_xlabel(r'$\Delta / \omega_m$', fontsize=20)
axes.set_ylabel(r'$g / \omega_m \times 10^{-5}$', fontsize=20)
# axes.set_zlabel(r'$E_N$', fontsize=14)
X, Y = np.meshgrid(Deltas / omega_m, g_0s / omega_m / 1e-5)
Z = np.array(E_N)
colors = ['b', 'w', 'r']
cmap = sns.diverging_palette(250, 15, s=75, l=40, n=11, center='light', as_cmap=True)
_sm = plt.cm.ScalarMappable(cmap=cmap, norm=mc.Normalize(vmin=min(_mini), vmax=max(_maxi)))
_sm.set_array([])
entan = axes.contourf(X, Y, Z, levels=11, cmap=cmap)
insta = axes.contourf(X, Y, np.array(Ins), cmap='gray')
# entan = axes.plot_wireframe(X, Y, Z, rstride=1, cstride=1, color='r', linewidth=0.25, antialiased=True)
# entan = axes.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
fig.colorbar(_sm, shrink=1, aspect=15)
# axes.view_init(30, 135)
plt.tight_layout()
plt.show()