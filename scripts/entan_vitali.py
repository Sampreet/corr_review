import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import scipy.linalg as sl

h_bar   = sc.hbar
K_B     = sc.k

w_b     = 2 * np.pi * 1e7
gamma   = 2 * np.pi * 1e2
lamb    = 810 * 1e-9
w_l     = 2 * np.pi * sc.c / lamb
w_c     = w_l + w_b
L       = 1e-3
finesse = 1.07e4
kappa   = 0.3 * w_b
m       = 5e-12
g_0     = w_c / L * np.sqrt(h_bar / m / w_b)
P       = 50e-3
print(kappa / w_b)
print(g_0)

Deltas = []
E = []
N = np.linspace(0, 5, 501)
for n in np.nditer(N):

    delta   = n * w_b
    wl      = w_c - delta
    T       = 40
    n_th    = 1 / (np.exp(h_bar * w_b / K_B / T - 1))
    E_0     = np.sqrt(2 * P * kappa / h_bar / wl)
    print(E_0)
    alpha   = np.abs(E_0 / (1j * delta + kappa))
    G       = g_0 * np.sqrt(2) * alpha
    A       = [ [0,     w_b,    0,      0       ],
                [-w_b,  -gamma, G,      0       ],
                [0,     0,      -kappa, delta   ],
                [G,     0,      -delta, -kappa  ]   ]
    D       = [ [0,     0,                  0,      0       ],
                [0,     gamma*(2*n_th+1),   0,      0       ],
                [0,     0,                  kappa,  0       ],
                [0,     0,                  0,      kappa   ]   ]
    V       =   sl.solve_lyapunov(np.array(A), -np.array(D))
    A       = [ [   V[0][0],    V[0][1]     ],
                [   V[1][0],    V[1][1]     ]   ]
    B       = [ [   V[2][2],    V[2][3]     ],
                [   V[3][2],    V[3][3]     ]   ]
    C       = [ [   V[0][2],    V[0][3]     ],
                [   V[1][2],    V[1][3]     ]   ]
    Sigma   = np.linalg.det(A) + np.linalg.det(B) - 2*np.linalg.det(C)
    eta_m   = np.sqrt((Sigma - np.sqrt(Sigma**2 - 4*np.linalg.det(V)))/2)
    en      = np.max([0, - np.log(2*eta_m)])

    Deltas.append(delta)
    E.append(en)
    
plt.figure()  
plt.plot(Deltas, E)
plt.show()