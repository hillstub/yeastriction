import dash
from dash import html
import dash_bootstrap_components as dbc

dash.register_page(__name__, path='/paper')

layout = dbc.Container([
    dbc.Row([
        dbc.Col(html.H2("Paper"), className="text-center my-4")
    ]),
    dbc.Row(
        dbc.Col([
            'Please cite',
            html.Ul([
                html.Li([
                    'Mans R, van Rossum HM, Wijsman M, Backx A, Kuijpers NGA, van den Broek M, Daran-Lapujade P, Pronk JT, van Maris AJA, Daran J-MG (2015) CRISPR/Cas9: a molecular Swiss army knife for simultaneous introduction of multiple genetic modifications in ',
                    html.I('Saccharomyces cerevisiae. '),
                    html.I('FEMS Yeast Research '),
                    html.B('15'),
                    '. ',
                    html.A('[link]', href="http://www.ncbi.nlm.nih.gov/pubmed/25743786", target="_blank")
                ])
            ])
        ], width=12)
    )
], fluid=True)
