import matplotlib.pyplot as plt
import numpy as np

omega_m = 1
# Delta_0 / omega_m
X       = np.linspace(-2, 2, 1001)
Kappa   = np.array([0.1, 0.2, 0.3, 0.4]) * omega_m
colors_1= ['#ed715f', '#d94732', '#a33424', '#6e291f']
colors_2= ['#337cf5', '#4a77c2', '#385b94', '#23324a']
plots   = []
legends = []

f = plt.figure()
plt.rcParams.update({'font.size': 12})
plt.xticks(np.arange(-3, 3, 0.5), fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel(r'$\Delta / \omega_m$', fontsize=14)

for i in range(len(Kappa)):
    kappa   = Kappa[i]
    Gamma   = []
    Shift   = []
    for j in range(len(X)):
        deno_1 = (kappa / omega_m)**2 / 4 + (1 - X[j])**2
        deno_2 = (kappa / omega_m)**2 / 4 + (1 + X[j])**2
        Gamma.append((kappa / omega_m) * (1 / deno_1 - 1 / deno_2))
        Shift.append(4 * ((1 - X[j]) / deno_1 - (1 + X[j]) / deno_2))
                
    plots.append(plt.plot(X, Gamma, c=colors_1[i]))
    plots.append(plt.plot(X, Shift, c=colors_2[i]))
    if i == len(Kappa) - 1:
        legends.append(r'$\gamma_{om} (\omega_m / g^2)$')
        legends.append(r'$4 \delta \omega_m  (\omega_m / g^2)$')

# plot
plt.legend(legends, loc='best')
plt.show()
