#!/usr/bin/env python

from pathlib import Path
import scanpy as sc
import scanpy_scripts as ss

# Arguments:
# 
# annData input
# bundle output dir



## Write matrices
# Call MTX-writing code used in scanpy-scripts
# Write .raw.X, .X and anything in layers, or as specified in .uns['atlas_bundle']['matrices']
# Add lines to MANIFEST for each matrix


## Write dimension reductions
# Write dimreds in same way as scanpy-scripts.
# Write everythin in .obsm, or as specified in .uns['atlas_bundle']['dimension_reductions']

# Write cell metadata 
# Write everything in 



# Have a 'prepare_for_atlas' function to take an annData object and record
# which slots will be exported to Atlas.

def prepare_for_atlas(
    x='filtered_normalised', 
    raw.x='raw', 
    layers='filtered'):


# Note which slots should be exported 

def tag_for_atlas(
    adata=None,
    type='matrix',
    slot=None,
    name=None,
    use_raw=False,
    file=None):

    adata.uns['scxa'][type][name] = {
        'slot' = slot,
        'use_raw' = use_raw,
        'file' = file
    }

def export_for_atlas(adata, dir = None):

    if dir is None:
        raise ValueError('No output directory provided')
    else:
        Path(dir).mkdir(parents=True, exist_ok=True)

    if not 'scxa' in adata.uns:
        raise ValueError('No manifest present for single-cell expression atlas export')

    for slot_type in adata.uns['scxa']:
        for name, manifest_entry in adata.uns[slot_type].items():
        
            if slot_type = 'matrix':
                if manifest_entry['slot'] != 'X':
                    layer = manifest_entry['slot']

                ss.cmd_utils.write_mtx(adata, fname_prefix = '%s/%s/' % (dir, manifest_entry['file']), use_raw=manifest_entry['use_raw'], layer=layer)

            if slot_type = 'obsm':
                ss.obj_utils.write_embedding(adata, key=slot, embed_fn = '%s/%s' % (dir, manifest_entry['file']))
                
            if slot_type = 'obs':
                ad.obs[['age', 'organism']]to_csv(embed_fn, sep=s, header=False, index=True)

            if slot_type = 'var':

    



    
    
