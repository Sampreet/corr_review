from matplotlib import colors as mc, cm
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import scipy.integrate as si
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
g_0     = 3.5e-5 * omega_m
E_0     = 6.4e3 * omega_m
E_1     = 3.2e3 * omega_m

Omega   = 2 * omega_m
tau     = 2 * np.pi / Omega
teval   = np.linspace(0, 50 * tau, 20001)

t = 0
v = [0, 0, 0]
c = [omega_m, Delta_0, gamma, kappa, g_0, E_0, E_1, Omega]
integrator = si.ode(model)
integrator.set_integrator('zvode')
integrator.set_initial_value(v, t)
integrator.set_f_params(c)

V2R = []
V2I = []
V0R = []
V1R = []
for i in range(1, len(teval)):
    t = teval[i]
    v = integrator.integrate(t)
    V2R.append(np.real(v[2]))
    V2I.append(np.imag(v[2]))
    V0R.append(np.real(v[0]))
    V1R.append(np.real(v[1]))

colors = sns.diverging_palette(250, 15, s=75, l=40, n=11, center='light', as_cmap=False)

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams.update({
    'font.family': 'Times New Roman',
    'font.style': 'normal',
    'font.variant': 'normal',
    'font.weight': 500,
    'font.stretch': 500,
    'font.size': 28,
    'axes.linewidth': 2
})

plt.figure()
plt.plot(V2R[:], V2I[:], color=colors[0], linewidth=1)
plt.plot(V2R[-1000:], V2I[-1000:], color=colors[-1], linewidth=5)
plt.xlabel(r'$Im[\langle a \rangle]$', labelpad=12, fontsize=32)
plt.ylabel(r'$Re[\langle a \rangle]$', labelpad=12, fontsize=32)
plt.xticks(fontsize=28)
plt.yticks(fontsize=28)
plt.ylim(-9000, -3000)

plt.tight_layout()

plt.figure()
plt.plot(V0R[:], V1R[:], color=colors[0], linewidth=1)
plt.plot(V0R[-1000:], V1R[-1000:], color=colors[-1], linewidth=5)
plt.xlabel(r'$\langle q \rangle$', labelpad=12, fontsize=32)
plt.ylabel(r'$\langle p \rangle$', labelpad=-24, fontsize=32)
plt.xticks(fontsize=28)
plt.yticks(fontsize=28)
plt.xlim(1750, 3500)
plt.ylim(-1500, 1500)

plt.tight_layout()

plt.show()

# from qom.ui import figure
# from qom.utils import axis

# _dim = 1001

# plot_params = {
#     'type': 'line',
#     'x_ticks': [1, 2, 5],
#     'show_legend': False,
#     'label_font_size': 28,
#     'tick_font_size': 24
# }

# x_dict = {
#     'var': 'xs',
#     'label': '$Im[\\langle a \\rangle]$',
#     'values': V2R
# }
# X = axis.StaticAxis(x_dict)

# y_dict = {
#     'var': 'ys',
#     'label': '$Re[\\langle a \\rangle]$',
#     'values': V2I
# }
# Y = axis.StaticAxis(y_dict)

# plotter = figure.Plotter(plot_params, X)

# X_v = axis.DynamicAxis([len(X.values)])
# X_v.values = X.values

# V = axis.DynamicAxis([len(Y.values)])
# V.values = Y.values

# plotter.update(X=X_v, Y=V)