from matplotlib import colors as mc
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import scipy.linalg as sl

omega_m = 2*np.pi*1e7                   
kappa   = 0.30 * omega_m
gamma   = 0.01 * omega_m
Delta   = - omega_m
G       = np.linspace(0, 0.03, 101) * omega_m
T       = np.linspace(0, 0.001, 101)

EN      = []
for g in np.nditer(G):
    temp = []
    for t in np.nditer(T):
        expr_1 = Delta**4*gamma*kappa + Delta**2*gamma**3*kappa/2 + Delta**2*gamma**2*kappa**2 + Delta**2*gamma*kappa**3/2 - 2*Delta**2*gamma*kappa*omega_m**2 + 4*Delta*g**2*gamma**2*omega_m + 8*Delta*g**2*gamma*kappa*omega_m + 4*Delta*g**2*kappa**2*omega_m + gamma**5*kappa/16 + gamma**4*kappa**2/4 + 3*gamma**3*kappa**3/8 + gamma**3*kappa*omega_m**2/2 + gamma**2*kappa**4/4 + gamma**2*kappa**2*omega_m**2 + gamma*kappa**5/16 + gamma*kappa**3*omega_m**2/2 + gamma*kappa*omega_m**4

        expr_2 = 4*Delta**2*gamma**2 + 16*Delta**2*omega_m**2 - 64*Delta*g**2*omega_m + gamma**2*kappa**2 + 4*kappa**2*omega_m**2

        if (expr_1 > 0 and expr_2 > 0):
            A = [   [   -kappa/2,   0,          0,          g           ],
                    [   0,          -kappa/2,   g,          0           ],
                    [   0,          g,         -gamma/2,    0           ],
                    [   g,          0,          0,          -gamma/2    ]   ]
            n = 1 / (np.exp(sc.hbar*omega_m/sc.k/t) - 1)
            D = [   [   -kappa/2,   0,          0,                  0                   ],
                    [   0,          -kappa/2,   0,                  0                   ],
                    [   0,          0,          -gamma*(2*n + 1)/2, 0                   ],
                    [   0,          0,          0,                  -gamma*(2*n + 1)/2  ]   ]
            V = sl.solve_lyapunov(np.array(A), -1*np.array(D))
            A = [   [   V[0][0],    V[0][1]     ],
                    [   V[1][0],    V[1][1]     ]   ]
            B = [   [   V[2][2],    V[2][3]     ],
                    [   V[3][2],    V[3][3]     ]   ]
            C = [   [   V[0][2],    V[0][3]     ],
                    [   V[1][2],    V[1][3]     ]   ]
            Sigma   = np.linalg.det(A) + np.linalg.det(B) - 2*np.linalg.det(C)
            eta_m   = np.sqrt((Sigma - np.sqrt(Sigma**2 - 4*np.linalg.det(V)))/2)
            en      = np.max([0, - np.log(2*eta_m)])
            temp.append(en)
        else:
            temp.append(0)
    
    EN.append(temp)
plt.figure()
plt.rcParams.update({'font.size': 12})
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel(r'$T (mK)$', fontsize=14)
plt.ylabel(r'$100g/\omega_m$', fontsize=14)
colors = [(1, 1, 1), (0.9, 0.75, 0.2)]
cmap = mc.LinearSegmentedColormap.from_list('hello', colors, 101)
plt.contourf(T*1000, G/omega_m*100, EN, cmap=cmap, alpha=1)
plt.colorbar()
plt.show()