import os
import dash
from dash import dcc, html, Input, Output, State, MATCH, ALL, callback
import pandas as pd
import dash_bootstrap_components as dbc
from dash import dash_table
from typing import List, Dict, Any
import logging
from sqlalchemy.exc import SQLAlchemyError

from models import Strain, Locus, Target, CRISPRSystem, search_targets, engine, LocalSession
from config import GENOMES_DIR


dash.register_page(__name__, path='/')

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def fetch_crispr_system_options() -> List[Dict[str, Any]]:
    """
    Fetch all available CRISPR systems from the database.

    Returns:
        List[Dict[str, Any]]: A list of dictionaries containing CRISPR system options.
        Each dictionary has 'value' (system ID) and 'label' (system name) keys.
    """
    try:
        with LocalSession() as session:
            crispr_systems = session.query(CRISPRSystem).all()
            options = [{'value': crispr_system.id, 'label': crispr_system.name} for crispr_system in crispr_systems]
        return options
    except SQLAlchemyError as e:
        logger.error(f"Error fetching CRISPR system options: {str(e)}")
        return []

def fetch_dna_build_options(crispr_system_id: int) -> List[Dict[str, str]]:
    """
    Fetch DNA build options for a specific CRISPR system.

    Args:
        crispr_system_id (int): The ID of the CRISPR system.

    Returns:
        List[Dict[str, str]]: A list of dictionaries containing DNA build options.
        Each dictionary has 'value' and 'label' keys representing the build method name.
    """
    try:
        with LocalSession() as session:
            crispr_system = session.query(CRISPRSystem).filter(CRISPRSystem.id == crispr_system_id).first()
            options = [{'value': oligo_build_method['name'], 'label': oligo_build_method['name']} for oligo_build_method in crispr_system.oligo_build_methods]
        return options
    except SQLAlchemyError as e:
        logger.error(f"Error fetching DNA build options: {str(e)}")
        return []

def fetch_strain_options() -> List[Dict[str, Any]]:
    """
    Fetch all available strains from the database.

    Returns:
        List[Dict[str, Any]]: A list of dictionaries containing strain options.
        Each dictionary has 'value' (strain ID) and 'label' (strain name) keys.
    """
    try:
        with LocalSession() as session:
            strains = session.query(Strain).order_by(Strain.name).all()
            options = [{'value': strain.id, 'label': strain.name} for strain in strains]
        return options
    except SQLAlchemyError as e:
        logger.error(f"Error fetching strain options: {str(e)}")
        return []

def fetch_locus_options(strain_id: int) -> List[Dict[str, Any]]:
    """
    Fetch locus options for a specific strain.

    Args:
        strain_id (int): The ID of the strain.

    Returns:
        List[Dict[str, Any]]: A list of dictionaries containing locus options.
        Each dictionary has 'value' (locus ID) and 'label' (locus ORF and symbol) keys.
    """
    try:
        with LocalSession() as session:
            loci = session.query(Locus).filter(Locus.strain_id == strain_id).order_by(Locus.orf).all()
            options = [{'value': locus.id, 'label': f'{locus.orf} ({locus.symbol})' if locus.symbol else locus.orf} for locus in loci]
        return options
    except SQLAlchemyError as e:
        logger.error(f"Error fetching locus options: {str(e)}")
        return []

layout = dbc.Container([
    html.H2('Yeastriction: Design guide RNAs for CRISPR-Cas9 in yeast', className='text-center my-4'),
    
    dbc.Card(
        [
            dbc.CardHeader('Inputs'),
            dbc.CardBody([
                dbc.Row([
                    dbc.Col([
                        dcc.Dropdown(
                            id='crispr-system-dropdown', 
                            options=fetch_crispr_system_options(), 
                            placeholder='Select a CRISPR system...',
                            className='mb-3'
                        ),
                    ], width=12)
                ]),
                dbc.Row([
                    dbc.Col([
                        dcc.Dropdown(
                            id='dna-build-method-dropdown', 
                            placeholder='Select a DNA build method...', 
                            className='mb-3'
                        ),
                    ], width=12)
                ]),
                dbc.Row([
                    dbc.Col([
                        dcc.Dropdown(
                            id='strain-dropdown', 
                            options=fetch_strain_options(), 
                            placeholder='Select a strain...',
                            className='mb-3'
                        ),
                    ], width=12)
                ]),
                dbc.Row([
                    dbc.Col([
                        dcc.Dropdown(
                            id='locus-dropdown', 
                            multi=True, 
                            placeholder='Select loci...', 
                            className='mb-3'
                        ),
                    ], width=12)
                ]),
            ])
        ],
        className='mb-4'
    ),
    
    dbc.Card([
            dbc.CardHeader('Loci'),
            dbc.CardBody([
                html.Div(id='targets-table-container')
            ]),
        ],
        className='mb-4'
    ),
    
    dbc.Card([
        dbc.CardHeader(['Selected Targets ', dcc.Clipboard(id="clipboard", style={"display": "inline-block"})]),
        dbc.CardBody([
            dash_table.DataTable(
                id='selected-targets-table',
                columns=[
                    {'name': 'Primer name', 'id': 'primer_name'},
                    {'name': 'Sequence', 'id': 'primer_sequence'}
                ],
                data=[],
                style_table={'overflowX': 'auto', 'boxShadow': '0 0 10px rgba(0, 0, 0, 0.1)', 'borderRadius': '5px'},
                style_header={
                    'backgroundColor': '#007bff', 
                    'color': 'white', 
                    'fontWeight': 'bold', 
                    'textAlign': 'left', 
                    'padding': '10px', 
                    'fontSize': '16px',
                    'borderTopLeftRadius': '5px',
                    'borderTopRightRadius': '5px'
                },
                style_cell={
                    'textAlign': 'left', 
                    'padding': '10px', 
                    'fontFamily': 'Arial', 
                    'fontSize': '14px',
                    'border': '1px solid #ddd',
                    'userSelect': 'text'  # Allows text selection
                },
                style_data_conditional=[
                    {
                        'if': {'row_index': 'odd'},
                        'backgroundColor': 'rgba(0, 123, 255, 0.1)'
                    },
                    {
                        'if': {'row_index': 'even'},
                        'backgroundColor': 'white'
                    },
                    {
                        'if': {'column_id': 'primer_sequence'},
                        'fontFamily': 'Courier New, monospace'
                    },                    
                ],
                style_as_list_view=True,
            )
        ])
        ],
        className='mb-4'
    ),
  
    dcc.Store(id='selected-targets-store', data={}),
], fluid=True)

@callback(
    Output('dna-build-method-dropdown', 'options'),
    Input('crispr-system-dropdown', 'value')
)
def set_dna_build_method_options(selected_crispr_system):
    """
    Callback to update DNA build method options based on the selected CRISPR system.

    Args:
        selected_crispr_system (int): The ID of the selected CRISPR system.

    Returns:
        List[Dict[str, str]]: A list of DNA build method options for the selected CRISPR system.
    """
    logger.info(f"Selected CRISPR system: {selected_crispr_system}")
    if selected_crispr_system is not None:
        return fetch_dna_build_options(selected_crispr_system)
    else:
        return []
    
@callback(
    Output('locus-dropdown', 'options'),
    Input('strain-dropdown', 'value')
)
def update_locus_dropdown(strain_id):
    """
    Callback to update locus dropdown options based on the selected strain.

    Args:
        strain_id (int): The ID of the selected strain.

    Returns:
        List[Dict[str, Any]]: A list of locus options for the selected strain.
    """
    logger.info(f"Selected strain: {strain_id}")
    if not strain_id:
        return []
    return fetch_locus_options(strain_id)

@callback(
    Output('targets-table-container', 'children'),
    Output('selected-targets-store', 'data', allow_duplicate=True),
    State('crispr-system-dropdown', 'value'),
    Input('locus-dropdown', 'value'),
    State('selected-targets-store', 'data'),
    prevent_initial_call=True
)
def find_targets(crispr_system_id, locus_ids, selected_targets):
    """
    Callback to find targets for selected loci and update the targets table.

    Args:
        crispr_system_id (int): The ID of the selected CRISPR system.
        locus_ids (List[int]): A list of selected locus IDs.
        selected_targets (Dict): Currently selected targets.

    Returns:
        Tuple[dash_bootstrap_components.Accordion, Dict]: 
            - An Accordion component containing target information for each locus.
            - Updated selected targets data.
    """

    if not locus_ids:
        return '', selected_targets
    accordions = []
    with LocalSession() as session:
        for locus_id in locus_ids:
            locus = session.query(Locus).filter(Locus.id == locus_id).first()
            targets = search_targets(session, locus_id, crispr_system_id=crispr_system_id)
            if not targets:
                accordions.append(dbc.AccordionItem(
                    title=f"{locus.display_name}",
                    children=html.Div('No targets found.'),
                    id={'type': 'accordion-item', 'index': locus_id}
                ))
            else:
                table_columns = [
                    {'name': 'Target sequence', 'id': 'sequence'},
                    {'name': 'AT content', 'id': 'at_content'},                    
                    {'name': 'RNA score', 'id': 'rna_fold_score'},
                    {'name': 'RNA structure', 'id': 'rna_fold_relevant_structure'},
                    {'name': 'z_score', 'id': 'z_score'},
                ]
                data = [{'id': t.id, 'orf': t.locus.orf, 'symbol': t.locus.symbol, 'strain_id': t.locus.strain_id, 'sequence': t.sequence, 'locus_id': t.locus_id, 'rna_fold_score': round(t.rna_fold.score,2), 'at_content': round(1-t.GC_content,2), 'z_score': round(t.z_score,2), 'rna_fold_relevant_structure': t.rna_fold.notation_binding_only} for t in targets]
                data = sorted(data, key=lambda x: x['z_score'], reverse=True)
                if str(locus_id) in selected_targets:
                    selected_row = next((index for (index, d) in enumerate(data) if d['id'] == selected_targets[str(locus_id)]['id']), 0)
                else:
                    selected_row = 0
                    selected_targets[str(locus_id)] = data[0]
                
                table = dash_table.DataTable(
                    columns=table_columns,
                    data=data,
                    row_selectable='single',
                    selected_rows=[selected_row],  # Pre-select the remembered row
                    id={'type': 'targets-table', 'index': locus_id},
                    style_table= {'overflowX': 'auto', 'maxHeight': '300px','overflowY': 'auto', 'boxShadow': '0 0 10px rgba(0, 0, 0, 0.1)', 'borderRadius': '5px'},
                    style_header={
                        'backgroundColor': '#007bff', 
                        'color': 'white', 
                        'fontWeight': 'bold', 
                        'textAlign': 'left', 
                        'padding': '10px', 
                        'fontSize': '16px',
                        'borderTopLeftRadius': '5px',
                        'borderTopRightRadius': '5px',
                        'position': 'sticky', 'top': '0','zIndex': 1
                    },
                    style_cell={
                        'textAlign': 'left', 
                        'padding': '10px', 
                        'fontFamily': 'Arial', 
                        'fontSize': '14px',
                        'border': '1px solid #ddd'
                    },
                    style_data_conditional=[
                        {
                            'if': {'row_index': 'odd'},
                            'backgroundColor': 'rgba(0, 123, 255, 0.1)'
                        },
                        {
                            'if': {'row_index': 'even'},
                            'backgroundColor': 'white'
                        },
                        {
                            'if': {'column_id': 'sequence'},
                            'fontFamily': 'Courier New, monospace'
                        },
                        {
                            'if': {'column_id': 'rna_fold_relevant_structure'},
                            'fontFamily': 'Courier New, monospace'
                        }                        
                    ],
                    style_as_list_view=True,
                )
                title = f"{locus.display_name} - {data[selected_row]['sequence']}"
                accordions.append(dbc.AccordionItem(
                    title=title,
                    children=table,
                    id={'type': 'accordion-item', 'index': locus_id},
                ))
    
    return dbc.Accordion(accordions, start_collapsed=True), selected_targets



@callback(
    Output({'type': 'accordion-item', 'index': MATCH}, 'title'),
    Input({'type': 'targets-table', 'index': MATCH}, 'selected_rows'),
    State({'type': 'targets-table', 'index': MATCH}, 'data'),
    State('selected-targets-store', 'data')
)
def update_selected_target(selected_rows, table_data, selected_targets):
    """
    Callback to update the title of an accordion item when a target is selected.

    Args:
        selected_rows (List[int]): Indices of selected rows in the targets table.
        table_data (List[Dict]): Data of the targets table.
        selected_targets (Dict): Currently selected targets.

    Returns:
        str: Updated title for the accordion item.
    """
    if selected_rows:
        with LocalSession() as session:
            selected_row = selected_rows[0]
            target_info = table_data[selected_row]
            locus_id = target_info['locus_id']
            selected_targets[str(locus_id)] = target_info
            locus = session.query(Locus).filter(Locus.id == locus_id).first()
            new_title = f"{locus.display_name} - {target_info['sequence']}"
            return new_title
    return dash.no_update

@callback(
    Output('selected-targets-store', 'data'),
    Input({'type': 'targets-table', 'index': ALL}, 'selected_rows'),
    State({'type': 'targets-table', 'index': ALL}, 'data'),
    State('selected-targets-store', 'data')
)
def store_selected_targets(selected_rows, table_data, selected_targets):
    """
    Callback to store selected targets in the selected-targets-store.

    Args:
        selected_rows (List[List[int]]): List of selected row indices for each targets table.
        table_data (List[List[Dict]]): Data of all targets tables.
        selected_targets (Dict): Currently selected targets.

    Returns:
        Dict: Updated selected targets data.
    """
    if selected_rows:
        for i, rows in enumerate(selected_rows):
            if rows:
                selected_row = rows[0]
                target_info = table_data[i][selected_row]
                locus_id = target_info['locus_id']
                selected_targets[str(locus_id)] = target_info
    return selected_targets


@callback(
    Output('selected-targets-table', 'data'),
    Input('dna-build-method-dropdown', 'value'),
    Input('selected-targets-store', 'data')
)
def update_selected_targets_table(dna_build_method, selected_targets):
    """
    Callback to update the selected targets table based on the DNA build method and selected targets.

    Args:
        dna_build_method (str): The selected DNA build method.
        selected_targets (Dict): Currently selected targets.

    Returns:
        List[Dict]: Updated data for the selected targets table.
    """
    rows = []
    with LocalSession() as session:
        
        for locus_id, target in selected_targets.items():
            locus = session.query(Locus).filter(Locus.id == locus_id).first()
            target = session.query(Target).filter(Target.id == target['id']).first()

            diagnostic_oligo_forward, diagnostic_oligo_reverse = locus.get_diagnostic_primers(session=session)
            
            for build_oligo in target.get_build_oligos(session=session, dna_build_method=dna_build_method):
                rows.append(build_oligo)

            repair_oligo_forward = {'primer_name': f'{locus.display_name}_repair oligo fw', 'primer_sequence': locus.repair_oligo_fw}
            repair_oligo_reverse = {'primer_name': f'{locus.display_name}_repair oligo rv', 'primer_sequence': locus.repair_oligo_rv}
            diagnostic_oligo_forward = {'primer_name': f'{locus.display_name}_dg fw', 'primer_sequence': diagnostic_oligo_forward}
            diagnostic_oligo_reverse = {'primer_name': f'{locus.display_name}_dg rv', 'primer_sequence': diagnostic_oligo_reverse}        
            rows.extend([
                repair_oligo_forward,
                repair_oligo_reverse,
                diagnostic_oligo_forward,
                diagnostic_oligo_reverse
            ])

    return rows

@callback(
    Output("clipboard", "content"),
    Input("clipboard", "n_clicks"),
    State('selected-targets-table', 'data'),
)
def selected(n, selected_targets):
    """
    Callback to copy selected targets data to clipboard.

    Args:
        n (int): Number of times the clipboard button has been clicked.
        selected_targets (List[Dict]): Data of the selected targets table.

    Returns:
        str: Tab-separated values of selected targets for clipboard.
    """
    df = pd.DataFrame(selected_targets)
    return df.to_csv(index=False, header=False, sep='\t')