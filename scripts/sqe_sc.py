from matplotlib import colors as mc
import matplotlib.pyplot as plt
import numpy as np

omega_m = 2*np.pi*1e7                         
kappa   = 2e-2 * omega_m
gamma   = 1e-5 * omega_m
E_0     = 4e11
Deltas  = np.linspace(-2, 2, 101) * omega_m
g_0s    = np.linspace(0, 2.5, 101) * 1e3

Bool  = []
for g_0 in np.nditer(g_0s):
    temp = []
    for Delta in np.nditer(Deltas):
        alpha = E_0 / np.sqrt(kappa**2 / 4 + Delta**2)
        g = g_0 * alpha

        expr_1 = Delta**4*gamma*kappa + Delta**2*gamma**3*kappa/2 + Delta**2*gamma**2*kappa**2 + Delta**2*gamma*kappa**3/2 - 2*Delta**2*gamma*kappa*omega_m**2 + 4*Delta*g**2*gamma**2*omega_m + 8*Delta*g**2*gamma*kappa*omega_m + 4*Delta*g**2*kappa**2*omega_m + gamma**5*kappa/16 + gamma**4*kappa**2/4 + 3*gamma**3*kappa**3/8 + gamma**3*kappa*omega_m**2/2 + gamma**2*kappa**4/4 + gamma**2*kappa**2*omega_m**2 + gamma*kappa**5/16 + gamma*kappa**3*omega_m**2/2 + gamma*kappa*omega_m**4

        expr_2 = 4*Delta**2*gamma**2 + 16*Delta**2*omega_m**2 - 64*Delta*g**2*omega_m + gamma**2*kappa**2 + 4*kappa**2*omega_m**2

        if (expr_1 > 0 and expr_2 > 0):
            temp.append(True)
        else:
            temp.append(False)
    
    Bool.append(temp)

plt.figure()
plt.rcParams.update({'font.size': 12})
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel(r'$\Delta / \omega_m$', fontsize=14)
plt.ylabel(r'$g (KHz)$', fontsize=14)
colors = [(1, 1, 1), (0.9, 0.75, 0.2)]
cmap = mc.LinearSegmentedColormap.from_list('hello', colors, 2)
plt.contourf(Deltas / omega_m, g_0s / 1e3, Bool, cmap=cmap, alpha=0.6)
plt.show()