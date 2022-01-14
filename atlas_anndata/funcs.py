import pkg_resources
import jsonschema
from jsonschema import validate
import yaml
from yaml import load, dump, Loader, Dumper
import os
import sys
import re

schema_file = pkg_resources.resource_filename('atlas_anndata', 'config_schema.yaml')

def load_doc(filename):
    with open(filename, "r") as stream:
    	try:
        		return yaml.safe_load(stream)
    	except yaml.YAMLError as exc:
        		print(exc)

def validate_config(config):

    schema = load_doc(schema_file) 

    try:
        validate(instance=config, schema=schema)
    except jsonschema.exceptions.ValidationError as err:
        print(err, file=sys.stderr)        
        return False
    return True

def check_slot(adata, slot_type, slot_name):

    print(f"Checking for {slot_type} {slot_name}")

    if slot_type == 'matrices':
        if slot_name == 'X':
            return hasattr(adata, 'X')
        elif slot_name == 'raw.X':
            return hasattr(adata.raw, 'X')
        else:
            return slot_name in adata.layers
    
    elif slot_type == 'dimension_reductions':
        return slot_name in adata.obsm

    elif slot_type == 'cell_groups':
        return slot_name in adata.obs.columns

    else:
        print(f"{slot_type} slot type not recognised", file=sys.stderr)
        return False 

def validate_anndata_with_config(config_file, anndata_file):

    config = load_doc(config_file)
    
    # First validate the anndata descripton file against the YAML schema

    print(f"Validating {config_file} against {schema_file}")
    config_status = validate_config(config)
    
    if (config_status):
        print("Config YAML file successfully validated")
    else:
        print("FAILURE")
        sys.exit(1)

    # Now check that the things the YAML said about the annData file are true

    print(f"Now checking config against anndata file")
    import scanpy as sc
    adata = sc.read(anndata_file)

    # Check that data is present at the locations indicated

    for slot_type in [ 'matrices', 'cell_groups', 'dimension_reductions' ]:
        if slot_type in config:
            for slot_def in config[slot_type]:
                if not check_slot(adata, slot_type, slot_def['slot']): 
                    print(f"{slot_type} entry {slot_def['slot']} not present in anndata file {anndata_file}")
                    sys.exit(1)

    print(f"annData file successfully validated against config {config_file}")
    return (config, anndata)

def obs_markers(adata, obs):
    markers_slot = f"markers_{obs}"

    if markers_slot in adata.uns.keys():
        return markers_slot
    else:
        return False

def slot_kind_from_name(slot_type, slot_name):

    from atlas_anndata import MISSING_STRING

    kind = MISSING_STRING
    search_map = {}

    if slot_type == "dimension_reductions":
        kind = f"{kind}: 'pca', 'tsne' or 'umap'"
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

