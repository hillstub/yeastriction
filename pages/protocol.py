import dash
from dash import html
import dash_bootstrap_components as dbc

dash.register_page(__name__, path='/protocol')

layout = dbc.Container([
    dbc.Row([
        dbc.Col(html.H2("Protocol Page"), className="text-center my-4")
    ]),
    dbc.Row([
        dbc.Col(html.P("This page provides a brief overview of the protocol for CRISPR-Cas9 genome editing in yeast. This method can be used for both single and multiple genetic modifications."), width=12)
    ]),
    dbc.Row([
        dbc.Col(html.H3("Protocol Overview"), className="my-4")
    ]),
    dbc.Row([
        dbc.Col(html.P("Method 1: In vivo assembly of CRISPR-Cas9 gRNA plasmid using the pMEL series. This method is efficient for editing a single locus and has a simple workflow prior to yeast transformation."), width=12)
    ]),
    dbc.Row([
        dbc.Col(html.P("Method 2: In vitro assembly of single and double gRNA plasmids using the pROS. This method is recommended for introducing multiple genetic modifications simultaneously. The pROS method uses gRNA plasmids containing two gRNA coding sequences, facilitating restriction at two loci."), width=12)
    ]),
    dbc.Row([
        dbc.Col(html.H4("Steps"), className="my-4")
    ]),
    dbc.Row([
        dbc.Col(html.Ol([
            html.Li("Design of the guideRNA (gRNA) primers using the Yeastriction tool or manually."),
            html.Li("Construction of the gRNA insert fragment by annealing complementary primers or using PCR."),
            html.Li("Construction of the linearised gRNA plasmid backbone via PCR."),
            html.Li("Construction of the double-stranded repair fragment for deletions, mutations, or DNA integration."),
            html.Li("Assembly of the gRNA expression plasmid using the NEBuilder reaction (for pROS method)."),
            html.Li("Transformation of Saccharomyces cerevisiae cells and selection of successful transformants."),
            html.Li("Confirmation and storage of the constructed plasmid via colony PCR or restriction analysis."),
            html.Li("Plasmid removal and troubleshooting.")
        ]), width=12)
    ]),
    dbc.Row([
        dbc.Col(
                html.P(["For detailed instructions, please check the supplementary materials in Mans R, van Rossum HM, et al. (2015) ",
                html.A('[link]', href="http://www.ncbi.nlm.nih.gov/pubmed/25743786", target="_blank")])
    , width=12)
    ])
], fluid=True)
