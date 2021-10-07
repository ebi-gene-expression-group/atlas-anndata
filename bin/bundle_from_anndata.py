#!/usr/bin/env python

import click
from pathlib import Path
import scanpy as sc
import scanpy_scripts as ss
import pathlib
import pandas as pd
import gzip
import os
from collections import OrderedDict
from scanpy_scripts.click_utils import  CommaSeparatedText
import re

@click.command()
@click.argument('anndata_file', type=click.Path(exists=True))
@click.argument('bundle_dir', type=click.Path(dir_okay=True))
@click.option('--droplet', default=False, help='Specify that this is a droplet experiment. There must be an obs parameter that can be used to differentiate cells from different runs.')
@click.option('--run-obs', default='run', help='The column in obs that differentiates cells from different runs. Will be used to separate run-wise and cell-wise metdata.')
@click.option('--exp-desc', default=None, help='Provide an experiment description file. If unspecified, the script will check for a slot called "experiment" in the .uns slot of the annData object, and use that to create a starting version of an IDF file.', type=click.Path(exists=True))
@click.option('--raw-matrix-slot', default=None, help='Slot where the most raw and least filtered expression data are stored. This is usually in raw.X, other values will be interpreted as layers (in .raw if prefixed with "raw". annData requires that .raw.X and .X match in .obs, so even when stored in .raw this will have to be the matrix afer cell filtering, but should ideally not be gene-filtered.')
@click.option('--filtered-matrix-slot', default=None, help='Layer name or "X", specifying storage location for gene-filtered matrix before transformations such as log transform, scaling etc.')
@click.option('--normalised-matrix-slot', default=None, help='Layer name or "X", specifiying storage location for normalised matrix.')
@click.option('--final-matrix-slot', default=None, help='Layer name or "X", specifiying storage location for the final matrix after processing. Will usually be in .X.')
@click.option('--write-obsms/--no--write--obsms', default=True, is_flag=True, help='Write dimension reductions from .obsm to the bundle?')
@click.option('--obsms', type=CommaSeparatedText(), default='all', help='Either "all" or a comma-separated list of obsm slots to write.')
@click.option('--atlas-style', default='run', help='Assume the tight conventions from SCXA, e.g. on .obsm slot naming')

def create_bundle(anndata_file, bundle_dir, droplet=False, run_obs='run', exp_desc=None, raw_matrix_slot=None, filtered_matrix_slot=None, normalised_matrix_slot=None, final_matrix_slot=None, write_obsms=True, obsms='all', atlas_style=True):
    """Build a bundle directory compatible with Single Cell Expression Atlas (SCXA) build proceseses
   
    \b 
    anndata_file - A file of the annData hdf5 specification, with all necessary information for SCXA.
    bundle_dir   - A directory in which to create the bundle.
    """

    adata = sc.read(anndata_file)
    pathlib.Path(bundle_dir).mkdir(parents=True, exist_ok=True)    

    manifest=_read_file_manifest(bundle_dir)

    if raw_matrix_slot is not None:
        print("Storing raw matrix")
        manifest = _write_matrix_from_adata(manifest, adata, slot=raw_matrix_slot, bundle_dir=bundle_dir, subdir='raw')        
    
    if filtered_matrix_slot is not None:
        print("Storing filtered matrix")
        manifest = _write_matrix_from_adata(manifest, adata, slot=filtered_matrix_slot, bundle_dir=bundle_dir, subdir='raw_filtered')         
    
    if normalised_matrix_slot is not None:
        print("Storing normalised matrix")
        manifest = _write_matrix_from_adata(manifest, adata, slot=normalised_matrix_slot, bundle_dir=bundle_dir, subdir='filtered_normalised')         
    
    if final_matrix_slot is not None:
        print("Storing normalised matrix")
        manifest = _write_matrix_from_adata(manifest, adata, slot=final_matrix_slot, bundle_dir=bundle_dir, subdir='final')         
       
    if write_obsms:
        print("Writing obsms to file")
        obsms = adata.obsm.keys() if obsms[0] == 'all' else obsms

        for obsm in obsms:
            print(f"Writing {obsm}")
            manifest = _write_obsm_from_adata(manifest, adata, obsm_name=obsm, bundle_dir=bundle_dir, atlas_style= atlas_style)

    # Write the final file manifest

    _write_file_manifest(bundle_dir, manifest) 

# Make filename and parameterisation from obsm names

def _write_obsm_from_adata(manifest, adata, obsm_name, bundle_dir, atlas_style = True):

    obsm_paramstring = obsm_name.replace('X_', '')
    filename = f"{obsm_paramstring}.tsv"
    
    description='embeddings'
    for embedding_type in [ 'umap', 'pca', 'tsne' ]:
        if embedding_type in obsm_paramstring:
            description=f'{embedding_type}_{description}'
            obsm_paramstring=re.sub(r"_?%s_?" % embedding_type, '', obsm_paramstring)
            break 
  
    # If following Atlas conventions, we can assume .obsm slots like
    # X_umap_neighbors_n_neighbors_20, which after above processing should now
    # be like neighbors_n_neighbors_20

    if atlas_style:
        parameterisation = obsm_paramstring.split('_')[-1]
    else:
        parameterisation = obsm_paramstring
    
    ss.obj_utils.write_embedding(adata, key=obsm_name, embed_fn=f"{bundle_dir}/{filename}")

    manifest = _set_manifest_value(manifest, description, filename, parameterisation)
    
    return manifest


def _read_file_manifest(bundle_dir):
    manifest_file=f"{bundle_dir}/MANIFEST"
    manifest=OrderedDict()

    if os.path.isfile(manifest_file):
        with open(manifest_file) as fp:
            header=fp.readline()
            for line in fp:
                line_parts=line.rstrip().split("\t")
                manifest = _set_manifest_value(manifest, line_parts[0], line_parts[1], line_parts[2])
            
    return manifest

def _write_file_manifest(bundle_dir, manifest):
    manifest_file=f"{bundle_dir}/MANIFEST"

    with open(manifest_file, 'w') as fh:
        fh.write('Description\tFile\tParameterisation\n')

        for description, v in manifest.items():
            for parameterisation, filename in v.items():
                fh.write(f"{description}\t{filename}\t{parameterisation}\n")

def _set_manifest_value(manifest, description, filename, parameterisation):
    if description not in manifest:
        manifest[description] = OrderedDict()
    
    manifest[description][parameterisation] = filename

    return manifest

def _write_matrix_from_adata(manifest, adata, slot, bundle_dir, subdir):
    
    layer = None
    use_raw = False

    if 'raw.' in slot:
        use_raw=True
        slot = slot.replace('raw.', '')
    
    if slot != 'X':
        layer = slot

    pathlib.Path(bundle_dir, subdir).mkdir(parents=True, exist_ok=True)    
    #ss.cmd_utils.write_mtx(adata, fname_prefix = '%s/%s/' % (bundle_dir, filename), use_raw=use_raw, use_layer=layer)
    write_mtx(adata, fname_prefix = f"{bundle_dir}/{subdir}/", use_raw=use_raw, use_layer=layer)
    
    for filename in [ 'matrix.mtx', 'barcodes.tsv', 'genes.tsv' ]:
        subfile=f"{subdir}/{filename}" 
        filepath=f"{bundle_dir}/{subfile}"
        with open(filepath, 'rb') as f_in, gzip.open("%s.gz" % filepath, 'wb') as f_out:
            f_out.writelines(f_in)
        os.remove(filepath)
    
    manifest = _set_manifest_value(manifest, 'mtx_matrix_content', f"{subdir}/matrix.mtx.gz", subdir)
    manifest = _set_manifest_value(manifest, 'mtx_matrix_cols', f"{subdir}/barcodes.tsv.gz", subdir)
    manifest = _set_manifest_value(manifest, 'mtx_matrix_rows', f"{subdir}/genes.mtx.gz", subdir)

    return manifest

def write_mtx(adata, fname_prefix='', var=None, obs=None, use_raw=False, use_layer=None):
    """Export AnnData object to mtx formt
    * Parameters
        + adata : AnnData
        An AnnData object
        + fname_prefix : str
        Prefix of the exported files. If not empty and not ending with '/' or '_',
        a '_' will be appended. Full names will be <fname_prefix>matrix.mtx,
        <fname_prefix>genes.tsv, <fname_prefix>barcodes.tsv
        + var : list
        A list of column names to be exported to gene table
        + obs : list
        A list of column names to be exported to barcode/cell table
    """
    if fname_prefix and not (fname_prefix.endswith('/') or fname_prefix.endswith('_')):
        fname_prefix = fname_prefix + '_'
    if var is None:
        var = []
    if obs is None:
        obs = []

    import scipy.sparse as sp
    if use_raw:
        var_source = adata.raw.var
        mat = sp.coo_matrix(adata.raw.X)
    else:
        var_source = adata.var
        if use_layer is not None:
            mat=sp.coo_matrix(adata.layers[use_layer])
        else:
            mat = sp.coo_matrix(adata.X)

    obs = list(set(obs) & set(adata.obs.columns))
    var = list(set(var) & set(var_source.columns))
    
    n_obs, n_var = mat.shape
    n_entry = len(mat.data)
    header = '%%MatrixMarket matrix coordinate real general\n%\n{} {} {}\n'.format(
        n_var, n_obs, n_entry)
    df = pd.DataFrame({'col': mat.col + 1, 'row': mat.row + 1, 'data': mat.data})
    mtx_fname = fname_prefix + 'matrix.mtx'
    gene_fname = fname_prefix + 'genes.tsv'
    barcode_fname = fname_prefix + 'barcodes.tsv'
    with open(mtx_fname, 'a') as fh:
        fh.write(header)
        df.to_csv(fh, sep=' ', header=False, index=False)

    obs_df = adata.obs[obs].reset_index(level=0)
    obs_df.to_csv(barcode_fname, sep='\t', header=False, index=False)
    var_df = var_source[var].reset_index(level=0)
    if not var:
        var_df['gene'] = var_df['index']
    var_df.to_csv(gene_fname, sep='\t', header=False, index=False)

if __name__ == '__main__':
    create_bundle()
