import numpy as np
import plotly.express as px
import plotly.graph_objects as go

data = {
    'Year': [2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020],
    'Web of Science': [7, 16, 16, 29, 27, 34, 46, 56, 59, 78, 36],
    'Scopus': [6, 26, 15, 32, 23, 35, 30, 46, 45, 47, 27],
    'color_map': {
        'Web of Science': 'green',
        'Scopus': 'blue'
    }
}
labels = {
    'Year': 'Year',
    'value': 'Number of Publications',
    'variable': 'Sources'
}

# fig = ff.create_distplot([W, S], group_labels)

# fig = px.bar(data_frame, x='Year', y=['Web of Science', 'Scopus'], color_discrete_map='color_map', labels=labels, barmode='group', opacity=0.8, width=1024, height=384)

fig = go.Figure(
    data = [
        # go.Bar(
        #     name = 'Web of Science', 
        #     x = data['Year'], 
        #     y = data['Web of Science'], 
        #     marker=dict(color='green', opacity=0.8)
        # ),
        go.Bar(
            name = 'Scopus', 
            x = data['Year'], 
            y = data['Scopus'], 
            text = data['Scopus'],
            marker=dict(
                color='brown', 
                opacity=0.8
            )
        )
    ],
    layout = go.Layout(
        barmode = 'group',
        plot_bgcolor = '#fff',
        legend_title = 'Sources',
        xaxis = dict(
            title = 'Year',
            showgrid = True, 
            gridcolor = '#ddd',
            tickmode = 'array',
            tickvals = data['Year'],
            ticktext = data['Year']
        ),
        yaxis = dict(
            title = 'Number of Publications',
            showgrid = True, 
            gridcolor = '#ddd'
        ),
        font = dict(
            family = 'Times New Roman',
            size = 16
        ),
        showlegend = False,
        width = 768, 
        height = 384
    )
)
fig.show()
