import numpy as np
import plotly.graph_objects as go

# squeezing parameter
r = 0.5
alpha = 1 + 1j
V_v = [[1/2, 0], [0, 1/2]]
V_c = [[1/2, 0], [0, 1/2]]
V_s = [[np.exp(-2*r) / 2, 0], [0, np.exp(2*r) / 2]]

qs = np.linspace(-5, 5, 501)
ps = np.linspace(-5, 5, 501)

W = []

for q in np.nditer(qs):

    temp = []
    for p in np.nditer(ps):
        # vacuum state
        X_v = np.array([[q], [p]])
        # displaced vacuum state
        X_d = np.array([[q - np.real(alpha)], [p - np.imag(alpha)]])
        # vacuum wigner function
        W_v = np.exp(-0.5 * np.transpose(X_v).dot(np.linalg.pinv(V_v)).dot(X_v)) / (2 * np.pi * np.sqrt(np.linalg.det(V_v)))
        # coherent wigner function
        W_c = np.exp(-0.5 * np.transpose(X_d).dot(np.linalg.pinv(V_c)).dot(X_d)) / (2 * np.pi * np.sqrt(np.linalg.det(V_c)))
        # squeezed wigner function
        W_s = np.exp(-0.5 * np.transpose(X_v).dot(np.linalg.pinv(V_s)).dot(X_v)) / (2 * np.pi * np.sqrt(np.linalg.det(V_s)))
        
        temp.append(W_v[0][0])

    W.append(temp)

project_z   = lambda x, y, z: z
colorsurf_z = project_z(qs, ps, W)
colorscale  = [ 
    (0.0, 'blue'), 
    (0.5, 'white'),
    (1.0, 'red')
]

data = [
    go.Surface(
        x=qs,
        y=ps,
        z=W,
        showlegend=False,
        showscale=False,
        colorscale=colorscale
    ),
    go.Surface(
        x=qs, 
        y=ps,
        z=list(0.5 * np.ones(np.array(W).shape)),
        showlegend=False,
        showscale=False,
        colorscale=colorscale,
        surfacecolor=colorsurf_z
    )
]

layout = go.Layout(
    autosize = False,
    scene = dict(
        xaxis_title='q',
        yaxis_title='p',
        zaxis_title='W(q, p)',
        camera_eye= dict(
            x=2, 
            y=1.5, 
            z=1.5
        )
    ),
    width = 512,
    height = 512
)

fig = go.Figure(
    data = data,
    layout = layout
)
# fig.show()

from qom.ui import figure
from qom.utils import axis
plot_params = {
    'z_label': '$W (q, p)$',
    'type': 'surface_cz',
    'label_font_size': 24,
    'tick_font_size': 20,
    'show_cbar': False
}

x_dict = {
    'var': 'xs',
    'label': '$q$',
    'values': qs
}
X = axis.StaticAxis(x_dict)

y_dict = {
    'var': 'ys',
    'label': '$p$',
    'values': ps
}
Y = axis.StaticAxis(y_dict)

plotter = figure.Plotter(plot_params, X, Y)

Z = axis.DynamicAxis([len(ps), len(qs)])
Z.values = W

plotter.update(Z=Z)