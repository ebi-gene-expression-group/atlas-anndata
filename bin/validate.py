#!/usr/bin/env python

import jsonschema
from jsonschema import validate
import yaml
from yaml import load, dump, Loader, Dumper
import click
import os
import sys

script_dir = os.path.dirname(os.path.realpath(__file__))

def load_doc(filename):
    with open(filename, "r") as stream:
    	try:
        		return yaml.safe_load(stream)
    	except yaml.YAMLError as exc:
        		print(exc)

def validate_config(config):
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

schema_file = os.path.join(script_dir, '..', 'config_schema.yaml')
schema = load_doc(schema_file) 

@click.command()
@click.argument('config_file', type=click.Path(exists=True))
@click.argument('anndata_file', type=click.Path(exists=True))

def validate_anndata(config_file, anndata_file):

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
                    print("%s entry %s not present in anndata file %s" % (slot_type, slot_def['slot'], anndata_file), file=sys.stderr)
                    sys.exit(1)

    print(f"annData file successfully validated against config {config_file}")

if __name__ == '__main__':
    validate_anndata()
