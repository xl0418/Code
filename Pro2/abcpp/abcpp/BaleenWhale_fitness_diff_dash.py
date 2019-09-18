import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
dir_path = 'c:/Liang/Googlebox/Research/Project2/BaleenWhales/'
test_file = dir_path + 'test_diff.csv'

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

spe_label = ['spe%i' % i for i in range(1,16)]


input_data = pd.read_csv(test_file)
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

        html.Div([
            html.H4('Select generation'),
            dcc.Dropdown(
                id='Gen_dropdown',
                options=[{'label': i, 'value': i} for i in POS],
                value=POS_default

            ),

        ],
            style={'width': '48%', 'display': 'inline-block'}),

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
    dash.dependencies.Output(component_id='heatmap', component_property='figure'),
    [dash.dependencies.Input(component_id='Gen_dropdown', component_property='value')]
)
def update_graph(Gen_dropdown):
    heatmap_data = \
    input_data[(input_data['Generation'] == Gen_dropdown)][
        ['Spe1', 'Spe2', 'Diff']]
    heatmap_data = pd.merge(data, heatmap_data, on=['Spe1', 'Spe2'], how='outer').fillna(0)
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
            title='MAJORITY OF SOLD AT ' + str.upper(
                Gen_dropdown) + ' IS ON ' + str.upper(maxsale['Spec1'][0]) + ' ' + str(
                maxsale['Spec2'][0])
        )

    }


if __name__ == '__main__':
    app.run_server(debug=False)