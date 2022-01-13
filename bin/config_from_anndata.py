#!/usr/bin/env python

# Create a starting config based on an anndata file

import click
import scanpy as sc
import yaml
from yaml import load, dump, Loader, Dumper
import re

missing = 'FILL ME'
string = f"{missing} with a string"
boolean = f"{missing} with a boolean"

@click.command()
@click.argument('anndata_file', type=click.Path(exists=True))
@click.argument('config_file', type=click.Path())
@click.option('--atlas-style', default=False, is_flag=True, help='Assume the tight conventions from SCXA, e.g. on .obsm slot naming?')
@click.option('--exp-name', required=True, help='Specify an Expression Atlas identifier that will be used for this experiment. If not set, a placeholder value E-EXP-1234 will be used and can be edited in the bundle later.')
@click.option('--droplet', default=False, is_flag=True, help='Is this a droplet experiment?')

def write_config(anndata_file, config_file, atlas_style = False, exp_name = None, droplet = False):
    
    adata = sc.read(anndata_file)

    config = {
        'exp-name': exp_name,
        'droplet': droplet,
        'matrices': [],
        'cell_groups': [],
        'dimension_reductions': list()
    }

    # Describe the matrices

    matrix_slots = []

    if hasattr(adata, 'raw') and hasattr(adata.raw, 'X'):
        matrix_slots.append('raw.X')

    matrix_slots = matrix_slots + list(adata.layers.keys())
    matrix_slots.append('X')

    for ms in matrix_slots:
        matrix_entry = {
            'slot': ms,
            'measure': 'counts',
            'cell_filtered': True,
            'gene_filtered': False,
            'normalised': False,
            'log_transformed': False,
            'scaled': True,
            'parameters': {}
        }
        
        # Use assumptions we can rely on for Atlas

        if atlas_style:
            if ms in ['filtered', 'normalised', 'X' ]:
                matrix_entry['gene_filtered'] = True
            if ms in ['normalised', 'X' ]:
                matrix_entry['normalised'] = True
            if ms == 'X':
                matrix_entry['log_transformed'] = True
                if droplet:
                    matrix_entry['scaled'] = True

        config['matrices'].append(matrix_entry)

    # Describe cell-wise metadata columns

    for obs in adata.obs.columns:

        obs_entry = {
            'slot': obs,
            'kind': slot_kind_from_name('cell_groups', obs),
            'parameters': extract_parameterisation('cell_groups', obs, atlas_style)
        }
        
        markers_slot = obs_markers(adata, obs)
        if markers_slot:
            obs_entry['markers'] = True
            obs_entry['markers_slot'] = markers_slot
        else:
            obs_entry['markers_slot'] = False    

        config['cell_groups'].append(obs_entry)
   
    # Describe dimension reductions stored in .obsm
 
    for obsm in adata.obsm.keys():
        obsm_entry = {
            'slot': obsm,
            'kind': slot_kind_from_name('dimension_reductions', obsm),
            'parameters': extract_parameterisation('dimension_reductions', obsm, atlas_style)
        }

        config['dimension_reductions'].append(obsm_entry)

    with open(config_file, 'w') as file:
        documents = yaml.dump(config, file)
  
def obs_markers(adata, obs):
    markers_slot = f"markers_{obs}"

    if markers_slot in adata.uns.keys():
        return markers_slot
    else:
        return False

def slot_kind_from_name(slot_type, slot_name):

    kind = string
    search_map = {}

    if slot_type == "dimension_reductions":
        kind = f"{string}: 'pca', 'tsne' or 'umap'"
        search_map = {
            '.*pca': 'pca', 
            '.*umap': 'umap',
            '.*tsne': 'tsne'
        }

    elif slot_type == 'cell_groups':
        kind = 'curation'        
        search_map = {
            'leiden': 'analysis',
            'log1p': 'analysis',
            '^n_': 'analysis',
            'total_': 'analysis',
            'pct_': 'analysis',
            'mito': 'analysis',
            'louvain': 'analysis'
        }
    
    for kind_pattern in search_map.keys():
        if re.match(r'{}'.format(kind_pattern), slot_name.lower()):
            kind = search_map[kind_pattern]
            break            

    return kind

# For Atlas exps there's some semantics to exploit in some slot names

def extract_parameterisation(slot_type, slot_name, atlas_style = False):
    
    parameters = {}

    if atlas_style:

        if slot_type == 'dimension_reductions':
            m = re.search(r'X_(umap|tsne|pca)_(.*)_(.*)', slot_name.replace('umap_neighbors', 'umap'))
            if m:
                parameters[m.group(2)] = string_to_numeric(m.group(3)) 
        
        elif slot_type == 'cell_groups':
            m = re.search(r'(louvain|leiden)_(.*)_(.*)', slot_name)
            if m:
                parameters[m.group(2)] = string_to_numeric(m.group(3)) 

    return parameters

# Convert to int or float

def string_to_numeric(numberstring):
    if numberstring.isdigit():
        return int(numberstring)
    else:
        return float(numberstring)

if __name__ == '__main__':
    write_config()
                       
