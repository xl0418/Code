import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'
test_file = dir_path + 'test_diff_full_ratio.csv'

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

spe_label = ['spe%i' % i for i in range(1,16)]


input_data = pd.read_csv(test_file)
input_data = input_data.drop(input_data.columns[[0]], axis=1)
data = []
for i in spe_label:
    for j in spe_label:
        data.append([i, j])
data = pd.DataFrame(data, columns=['Spe1', 'Spe2'])

POS = input_data['Generation'].unique()
POS_default = input_data['Generation'][400000]
df_default = input_data[
    (input_data['Generation'] == POS_default)]

app.layout = html.Div([
    html.H1('Difference in fitness'),
    html.Div([

        # html.Div([
        #     html.H4('Select generation'),
        #     dcc.Dropdown(
        #         id='Gen_dropdown',
        #         options=[{'label': i, 'value': i} for i in POS],
        #         value=POS_default
        #
        #     ),
        #
        # ],
        #     style={'width': '48%', 'display': 'inline-block'}),

        html.Div([
            dcc.Slider(
                id='Gen_dropdown',
                min=0,
                max=POS.max(),
                step=1,
                value=10,
                marks={
                    0: {'label': 'Start', 'style': {'color': '#77b0b1'}},
                    4430: {'label': '1'},
                    21684: {'label': '2'},
                    24118: {'label': '3'},
                    25302: {'label': '4'},
                    26992: {'label': '5'},
                    29419: {'label': '6'},
                    29929: {'label': '7'},
                    31032: {'label': '8'},
                    35134: {'label': '9'},
                    35758: {'label': '10'},
                    37629: {'label': '11'},
                    39287: {'label': '12'},
                    40364: {'label': '13'},
                    41168: {'label': 'P', 'style': {'color': '#f50'}}
                }
            ),
            html.Div(id='slider-output-container', style={'margin-top': 20})
        ],
            style={'width': '80%', 'display': 'inline-block'}),


        dcc.Graph(id='heatmap',
                  figure={
                      'data': [go.Heatmap(
                          x=df_default['Spe1'],
                          y=df_default['Spe2'],
                          z=df_default['Diff'],
                          name='first legend group',
                          colorscale='Viridis')],
                      'layout': go.Layout(
                          xaxis=dict(title='Species'),
                          yaxis=dict(title='Species'),
                      )

                  })
    ]),

])

@app.callback(
    dash.dependencies.Output('slider-output-container', 'children'),
    [dash.dependencies.Input('Gen_dropdown', 'value')])
def update_output(value):
    return 'You have selected the time point {}'.format(value)

@app.callback(
    dash.dependencies.Output(component_id='heatmap', component_property='figure'),
    [dash.dependencies.Input(component_id='Gen_dropdown', component_property='value')]
)
def update_graph(Gen_dropdown):
    heatmap_data = \
    input_data[(input_data['Generation'] == Gen_dropdown)][
        ['Spe1', 'Spe2', 'Diff']]
    # heatmap_data = pd.merge(data, heatmap_data, on=['Spe1', 'Spe2'], how='outer').fillna(0)
    print(Gen_dropdown)
    maxsale = heatmap_data[heatmap_data['Diff'] == heatmap_data['Diff'].max()]
    maxsale = maxsale.reset_index()
    print(maxsale)
    return {
        'data': [go.Heatmap(
            x=heatmap_data['Spe1'],
            y=heatmap_data['Spe2'],
            z=heatmap_data['Diff'],
            xgap=2,
            ygap=2,
            colorscale='Viridis')],
        'layout': go.Layout(
            title=''
        )

    }


if __name__ == '__main__':
    app.run_server(debug=False)