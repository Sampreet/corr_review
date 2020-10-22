import numpy as np
import plotly.graph_objects as go
import scipy.constants as sc
import scipy.linalg as sl
import seaborn as sns

omega_m = 2*np.pi*1e7                   
kappa   = 7.5e-1 * omega_m
gamma   = 1e-5 * omega_m
T       = 0.4e-3
# g_0s    = np.array([2.5e3, 1.5e3])
g_0s    = np.array([4e-5, 3e-5]) * omega_m
opacities=[0.8, 1]
E_0     = 6.4e3 * omega_m
Delta   = omega_m
legends = [r'$g / \omega_m = 0.15$', r'$g / \omega_m = 0.25$']
colors  = sns.diverging_palette(250, 15, s=75, l=40, n=11, center='light', as_cmap=False)
print(str(colors[-1]))
colors  = ['rgb' + str(colors[-1]), 'rgb' + str(colors[0])]
X       = [r'$V_{11}$', r'$V_{22}$', r'$V_{33}$', r'$V_{44}$', r'$V_{12}$', r'$V_{34}$', r'$V_{13}$', r'$-V_{24}$', r'$V_{14}$', r'$V_{23}$']
data    = []

for i in range(len(g_0s)):
    g_0 = g_0s[i]
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
        print(V)

        data.append(
            go.Bar(
                name = legends[i],
                x = X,
                y = [V[0][0], V[1][1], V[2][2], V[3][3], V[0][1], V[2][3], V[0][2], -V[1][3], V[0][3], V[1][2]],
                marker = dict(
                    color = colors[i],
                    opacity = opacities[i]
                ),
                orientation = 'v'
            )
        )

fig = go.Figure(
    data = data,
    layout = go.Layout(
        barmode = 'overlay',
        plot_bgcolor = '#fff',
        xaxis = dict(
            title = '',
            showgrid = True, 
            gridcolor = '#fff',
            color = 'black'
        ),
        yaxis = dict(
            title = '',
            showgrid = True, 
            gridcolor = '#fff',
            tickvals = [0, 0.25, 0.5],
            tickfont = dict(
                size = 16
            ),
            color = 'black'
        ),
        font = dict(
            family = 'Times New Roman',
            size = 20
        ),
        showlegend = False,
        width = 512, 
        height = 256
    )
)
fig.show()