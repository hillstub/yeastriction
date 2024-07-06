import os
import sys
import time
from datetime import datetime
import base64
import io
import subprocess
import threading

from models import Strain, Locus, Target, CRISPRSystem, search_targets

from sqlalchemy import create_engine, select, delete, update
from sqlalchemy.orm import sessionmaker

import dash
from dash import html
import dash_bootstrap_components as dbc
from dash import dcc, html, Input, Output, State, MATCH, ALL, callback

from pathlib import Path
import pandas as pd

DATABASE_URL = 'sqlite:///./data/database.db'

GENOMES_DIR = os.getenv('GENOMES_DIR')

engine = create_engine(DATABASE_URL)
Session = sessionmaker(bind=engine)

dash.register_page(__name__, path='/import')

if '--allow-import' in sys.argv:
    # rewrite the function to use pandas
    def import_loci(strain_name: str, df: pd.DataFrame):
        with Session() as session:
            strain = session.execute(select(Strain).where(Strain.name == strain_name)).scalars().first()
            
            if not strain:
                strain = Strain(name=strain_name)
                session.add(strain)
                session.commit()
                session.refresh(strain)
            
            # Remove loci associated with the strain
            session.execute(delete(Locus).where(Locus.strain_id == strain.id))
            session.commit()
            
            # Read the .tab file
            #file_path = Path(GENOMES_DIR) / f'{strain_name}.tab'
            #df = pd.read_csv(file_path, sep='\t', header=None, names=['orf', 'symbol', 'start_orf', 'end_orf', 'sequence'])
            
            # if there are multiple rows with the same orf, only take the first one
            df = df.drop_duplicates(subset=['orf'], keep='first')

            df['strain_id'] = strain.id
            df['start_orf'] = df['start_orf'] - 1
            df['symbol'] = df['symbol'].replace('', None)
            df['created'] = datetime.now()
            
            # create Locus objects from the dataframe
            loci = [Locus(**row.to_dict()) for _, row in df.iterrows()]
            session.add_all(loci)
            session.commit()

            # get non-empty non-unique symbols
            non_unique_symbols =  df['symbol'].dropna().loc[df['symbol'].duplicated()].tolist()
            if len(non_unique_symbols) > 0:
                print(non_unique_symbols)
                #session.update(Locus).where(Locus.strain_id == strain.id, Locus.symbol.in_(non_unique_symbols)).values(symbol=None)
                session.execute(update(Locus).where(Locus.strain_id == strain.id, Locus.symbol.in_(non_unique_symbols)).values(symbol=None))
                session.commit()
            
            return {"message": "Import complete"}
        
    def save_fasta_file(content, filename):
        data = content.split(',')[1]
        file_data = base64.b64decode(data)
        file_path = Path(GENOMES_DIR) / filename
        with open(file_path, 'wb') as f:
            f.write(file_data)
        return file_path

    def index_fasta_file(file_path):
        try:
            index_base = file_path.with_suffix('')  # Remove the .fasta or .fa extension
            subprocess.run(['bowtie-build', str(file_path), str(index_base)], check=True)
            return {"message": f"Indexing of {file_path} complete"}
        except subprocess.CalledProcessError as e:
            return {"message": f"Error indexing {file_path}: {e}"}

    def parse_contents(contents, filename):
        content_type, content_string = contents.split(',')
        try:
            if filename.endswith('.tab'):
                decoded = base64.b64decode(content_string)
                df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep='\t', header=None, names=['orf', 'symbol', 'start_orf', 'end_orf', 'sequence'])
                strain_name = Path(filename).stem
                result = import_loci(strain_name, df)
                return html.Div([f'File {filename}: {result["message"]}'])
            elif filename.endswith('.fasta') or filename.endswith('.fa'):
                file_path = save_fasta_file(contents, filename)
                result = index_fasta_file(file_path)
                return html.Div([f'File {filename}: {result["message"]}'])
            else:
                return html.Div(['Unsupported file format'])
        except Exception as e:
            return html.Div([f'There was an error processing the file {filename}: {str(e)}'])

    layout = dbc.Container([
        html.H2('Genome Importer', className='text-center my-4'),
        dbc.Row([
            dbc.Col(dcc.Upload(
                id='upload-data',
                children=html.Div(['Drag and Drop or ', html.A('Select Files')]),
                style={
                    'width': '100%',
                    'height': '60px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '10px'
                },
                multiple=True
            ), width=12),
        ]),
        dbc.Row([
            dbc.Col(dbc.Progress(id='progress-bar', value=0, max=10, style={'width': '100%', 'visibility': 'hidden'}), width=12)
        ]),
        dbc.Row([
            dbc.Col(html.Div(id='output-data-upload'), width=12)
        ])
    ], fluid=True)

    @callback(
        output=Output('output-data-upload', 'children'),
        inputs=[Input('upload-data', 'contents'), State('upload-data', 'filename')],
        running=[
            (Output('upload-data', 'disabled'), True, False),
            (Output('progress-bar', 'style'), {'visibility': 'visible'}, {'visibility': 'hidden'})
        ],
        progress=[Output('progress-bar', 'value'), Output('progress-bar', 'max'), Output('output-data-upload', 'children')],
        prevent_initial_call=True,
        background=True
    )
    def update_output(set_progress, contents, filenames):
        if contents is not None:
            file_dict = {}
            for content, filename in zip(contents, filenames):
                strain_name = Path(filename).stem
                ext = Path(filename).suffix
                if strain_name not in file_dict:
                    file_dict[strain_name] = {}
                file_dict[strain_name][ext] = content
            
            children = []
            total_files = len(filenames)
            processed_files = 0

            for strain_name, files in file_dict.items():


                if '.tab' in files and ('.fasta' in files or '.fa' in files):
                    fasta_key = '.fasta' if '.fasta' in files else '.fa'
                    tab_content = files['.tab']
                    fasta_content = files[fasta_key]

                    children.append(parse_contents(tab_content, f'{strain_name}.tab'))
                    processed_files += 1
                    set_progress((processed_files, total_files, children))                    

                    children.append(parse_contents(fasta_content, f'{strain_name}{fasta_key}'))
                    processed_files += 1
                    set_progress((processed_files, total_files, children))                    
                else:
                    processed_files += 2
                    children.append(html.Div([f'Error: Missing .tab or .fasta file for strain {strain_name}']))
                    set_progress((processed_files, total_files, children)) 


            return children
else:
    layout = dbc.Container([
        dbc.Row([
            dbc.Col(html.H2("Access Denied"), className="text-center my-4")
        ]),
        dbc.Row([
            dbc.Col(html.P("You do not have permission to access this page. Please run the script with the --allow-import flag."), className="mb-4")
        ])
    ], fluid=True)