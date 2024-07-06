import os
import sys
import dash
import dash_bootstrap_components as dbc
from dash import html, dcc, DiskcacheManager
from dotenv import load_dotenv
import diskcache
from models import initialize_database

load_dotenv()

cache = diskcache.Cache("./cache")
background_callback_manager = DiskcacheManager(cache)

# Create the database if it doesn't exit yet
initialize_database()

# Initialize the Dash app
app = dash.Dash(__name__, use_pages=True, external_stylesheets=[dbc.themes.BOOTSTRAP], suppress_callback_exceptions=True, background_callback_manager=background_callback_manager)

# Define the app layout
app.layout = html.Div(children=[
    dcc.Location(id='url', refresh=False),
    dbc.NavbarSimple(
        brand='Yeastriction',
        brand_href='/',
        color='primary',
        dark=True,
        sticky='top',
        children=[
            dbc.NavItem(dbc.NavLink('Import', href='/import')) if '--allow-import'in sys.argv else None,            
            dbc.NavItem(dbc.NavLink('Protocol', href='/protocol')),
            dbc.NavItem(dbc.NavLink('Paper', href='/paper')),
            dbc.NavItem(dbc.NavLink('Github', href='https://github.com/your-repo', target='_blank'))
        ]
    ),
    dash.page_container  # This will contain the content of the pages
])

if __name__ == '__main__':
    app.run_server(host='0.0.0.0', debug=True)
