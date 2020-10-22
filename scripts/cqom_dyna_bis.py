import matplotlib.pyplot as plt
import numpy as np

omega_m = 1                             
kappa   = 0.30 * omega_m     
gamma   = 0.01 * omega_m
A_l     = 50 * kappa

# Delta_0 / omega_m
arr_x   = np.linspace(-1.5, 1.5, 1001)
G_0     = np.array([0.0, 0.001, 0.002, 0.003, 0.005, 0.008]) * omega_m
colors  = ['#402f2d', '#542f29', '#6e291f', '#a33424', '#d94732', '#ed715f']
X       = []
N       = []
plots   = []
legends = []

f = plt.figure()
plt.rcParams.update({'font.size': 12})
plt.xticks(np.arange(-2, 2.5, 0.5), fontsize=12)
plt.yticks(fontsize=12)
plt.ylabel(r'$N / N_{max}$', fontsize=14)
plt.xlabel(r'$\Delta_0 / \omega_m$', fontsize=14)

for i in range(len(G_0)):
    X = []
    N = []
    g_0 = G_0[i]
    index_start = 0
    index_stop = 0
    for j in range(len(arr_x)):
        x = arr_x[j]
        a = g_0**2 / omega_m / kappa
        b = x * omega_m / kappa
        c = gamma / omega_m
        d = A_l / kappa
        temp = (c**2 / 4  + 1)
        coeffs = []
        coeffs.append(4 * a**2)
        coeffs.append(- 4 * a * b * temp)
        coeffs.append((b**2 + 1 / 4) * temp**2)
        coeffs.append(- d**2 * temp**2)
        roots = np.roots(coeffs)

        count = 0
        for root in np.nditer(roots):
            if np.imag(root) == 0.0:
                X.append(x)
                N.append(np.real(root))
                count += 1
        if count > 1:
            if index_start == 0:
                index_start = j
                index_stop = j
            else:
                index_stop = j
                
    N = N/max(N)
    if index_start != 0:
        plt.axvspan(arr_x[index_start], arr_x[index_stop], color=colors[i], alpha=0.5)
    plots.append(plt.scatter(X, N, s=5, c=colors[i], marker='o', alpha=0.75))
    legends.append(r'$g_0 / \omega_m = ' + str(g_0) + '$')

# plot
plt.legend(plots, legends, loc='best')
plt.show()
