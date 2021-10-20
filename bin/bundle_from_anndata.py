#!/usr/bin/env python

import click
from pathlib import Path
import shutil
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
@click.option('--write-obsms/--no-write-obsms', default=True, is_flag=True, help='Write dimension reductions from .obsm to the bundle?')
@click.option('--obsms', type=CommaSeparatedText(), help='A comma-separated list of obsm slots to write. Default is to write them all.')
@click.option('--write-clusters/--no--write--clusters', default=True, is_flag=True, help='Write cluster data to the bundle?')
@click.option('--clusters', type=CommaSeparatedText(), help='A comma-separated list of .obs slots corresponding to unsupervised clusterings.')
@click.option('--clusters-field-pattern', default='louvain', help='If --clusters not supplied, a string to be used to select columns from .obs representing unsupervised clusterings.')
@click.option('--default-clustering', help='Where --write-clusters is set, which clustering should be set as the default? If not set, the middle (or first middle) clustering will be selected, or if --atlas-style is set, this will be the clustering corresponding to a resolution of 1.')
@click.option('--write-markers/--no--write--markers', default=True, is_flag=True, help='Write marker data to the bundle?')
@click.option('--marker-clusterings', type=CommaSeparatedText(), help='A comma-separated list of clusterings for which to write markers. marker results are expected to be stored in .uns under keys like "markers_{clustering}. Defaults to all selected clusterings.')
@click.option('--metadata-marker-fields', type=CommaSeparatedText(), help='A comma-separated list of .obs slots corresponding to metadata variables for which markers have been derived.')
@click.option('--atlas-style', default='run', help='Assume the tight conventions from SCXA, e.g. on .obsm slot naming')
@click.option('--gene-name-field', default='gene_name', help='Field in .var where gene name (symbol) is stored.')

def create_bundle(anndata_file, bundle_dir, droplet=False, run_obs='run', exp_desc=None, raw_matrix_slot=None, filtered_matrix_slot=None, normalised_matrix_slot=None, final_matrix_slot=None, write_obsms=True, obsms=None, write_clusters = True, clusters = None, clusters_field_pattern = 'louvain', default_clustering = None, write_markers = True, marker_clusterings=None, metadata_marker_fields=None, atlas_style=True, gene_name_field='gene_name'):
    """Build a bundle directory compatible with Single Cell Expression Atlas (SCXA) build proceseses
   
    \b 
    anndata_file - A file of the annData hdf5 specification, with all necessary information for SCXA.
    bundle_dir   - A directory in which to create the bundle.
    """

    adata = sc.read(anndata_file)

    # For any genes without names, assign the ID
    genes_with_missing_names = list(adata.var_names[pd.isnull(adata.var[gene_name_field])])
    adata.var[gene_name_field] = adata.var[gene_name_field].cat.add_categories(genes_with_missing_names)
    adata.var.loc[genes_with_missing_names, gene_name_field] = genes_with_missing_names

    if Path(bundle_dir).is_dir():
        shutil.rmtree(bundle_dir)
    
    pathlib.Path(bundle_dir).mkdir(parents=True)    

    manifest=_read_file_manifest(bundle_dir)

    if raw_matrix_slot is not None:
        print("Storing raw matrix")
        manifest = _write_matrix_from_adata(manifest, adata, slot=raw_matrix_slot, bundle_dir=bundle_dir, subdir='raw', gene_name_field=gene_name_field)        
    
    if filtered_matrix_slot is not None:
        print("Storing filtered matrix")
        manifest = _write_matrix_from_adata(manifest, adata, slot=filtered_matrix_slot, bundle_dir=bundle_dir, subdir='raw_filtered', gene_name_field=gene_name_field)         
    
    if normalised_matrix_slot is not None:
        print("Storing normalised matrix")
        manifest = _write_matrix_from_adata(manifest, adata, slot=normalised_matrix_slot, bundle_dir=bundle_dir, subdir='filtered_normalised', gene_name_field=gene_name_field)         
    
    if final_matrix_slot is not None:
        print("Storing final matrix")
        manifest = _write_matrix_from_adata(manifest, adata, slot=final_matrix_slot, bundle_dir=bundle_dir, subdir='final', gene_name_field=gene_name_field)         
       
    if write_obsms:
        print("Writing obsms to file")
        obsms = adata.obsm.keys() if obsms is None else obsms

        for obsm in obsms:
            print(f"Writing {obsm}")
            manifest = _write_obsm_from_adata(manifest, adata, obsm_name=obsm, bundle_dir=bundle_dir, atlas_style= atlas_style)
    
    if write_clusters:
        print("Writing clusterings to file")
 
        # If no list of cluster variables is specified, select automatically
        # using the field pattern. If the fields are supplied, use the ones
        # that are actually present in .obs.
       
        if clusters is None:
            clusters = [ x for x in adata.obs.columns if clusters_field_pattern in x ]
        else:
            clusters = list(set(clusters) & set(adata.obs.columns))
        
        if len(clusters) == 0:
            print(f"No .obs slots matching clusters specifications.") 
            return False

        manifest = _write_clusters_from_adata(manifest, adata, clusters = clusters, bundle_dir = bundle_dir, atlas_style = atlas_style)       
        
    if write_markers:
        print("Writing markers to file")
        manifest = _write_markers_from_adata(manifest, adata, clusters = clusters, marker_clusterings = marker_clusterings, metadata_marker_fields = metadata_marker_fields, bundle_dir = bundle_dir, atlas_style = atlas_style)       
 
    # Write the final file manifest

    _write_file_manifest(bundle_dir, manifest) 

# Write dimension reductions from .obsm slots

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

def _select_clusterings(adata, clusters, atlas_style = True):
    
    clusterings = list(dict(sorted(dict(zip( clusters, [ abs(1 - float(c.split('_')[-1])) for c in clusters ])).items(), key=lambda x: x[1])).keys())

    clustering_to_nclust = dict(zip( clusterings, [ len(adata.obs[c].unique()) for c in clusterings ]))

    if atlas_style:

        # Only keep the first marker set for a given k. For Atlas, ranked as above,
        # this will be the makers from the resolution closest to 1 of clashing
        # sets.

        kept_clusterings = {}
        for k, v in clustering_to_nclust.items():
            if v not in kept_clusterings.values():
                kept_clusterings[k] = v

        clustering_to_nclust = kept_clusterings

    return dict(sorted(clustering_to_nclust.items(), key=lambda item: item[1]))

def _write_clusters_from_adata(manifest, adata, clusters, bundle_dir, atlas_style = True, default_clustering = None):
   
    if default_clustering is None:
        if atlas_style:
            default_clustering = 'louvain_resolution_1.0'
        else:
            default_clustering = clusterings[math.floor((len(clusterings)-1)/2)]
     
    clustering_to_k = _select_clusterings(adata, clusters = clusters, atlas_style = atlas_style)

    # Write the complete clusters file in order
    clusters_output = open(f"{bundle_dir}/clusters_for_bundle.txt", mode="w")
    
    with open(f"{bundle_dir}/clusters_for_bundle.txt", mode="w") as fp:
        fp.write("sel.K\tK\t"+"\t".join(list(adata.obs_names))+"\n")

        for clustering, k in clustering_to_k.items():
            selected = 'TRUE' if clustering == default_clustering else 'FALSE'
            cluster_memberships = list(adata.obs[clustering])
            if min([ int(x) for x in cluster_memberships ]) == 0:
                cluster_memberships = [ str(int(x) + 1) for x in cluster_memberships ]

            fp.write("\t".join([selected, str(k)]))
            fp.write("\t"+"\t".join(list(cluster_memberships)))
            fp.write("\n")
            
    manifest = _set_manifest_value(manifest, 'cluster_memberships', f"{bundle_dir}/clusters_for_bundle.txt", '')

    return manifest

# For Atlas we store marker sets by integer number of clusters. Where cluster
# nubmers clash between marker sets (e.g. for different resolution values) we
# use the set from a resolution closest to 1
# TODO: fix this for cell type markers

def _write_markers_from_adata(manifest, adata, clusters, marker_clusterings = None, metadata_marker_fields = None, bundle_dir = None, atlas_style = True):
   
    clustering_to_k = _select_clusterings(adata, clusters = clusters, atlas_style = atlas_style)
   
    # If no explicit marker sets suppplied, select those corresponding to the
    # selected clusterings (where available)
 
    marker_groupings = []

    if marker_clusterings is None:
        marker_groupings = [ x for x in clustering_to_k.keys() if f"markers_{x}" in adata.uns.keys() ]
   
    if metadata_marker_fields is not None:
        marker_groupings = marker_groupings + metadata_marker_fields

    missing_marker_sets = [ x for x in marker_groupings if f"markers_{x}" not in adata.uns.keys() ] 
    if len(missing_marker_sets) > 0:
       raise Exception("Some supplied marker clusterings do not have marker results in .uns: %s" % ','.join(missing_marker_sets))         

    for mg in marker_groupings:
        marker = f"markers_{mg}"

        de_tbl = ss.lib._diffexp.extract_de_table(adata.uns[marker])
        de_tbl = de_tbl.loc[de_tbl.genes.astype(str) != 'nan', :]

        # Reset cluster numbering to be from 1 if required

        if de_tbl['cluster'].min() == '0':
            de_tbl['cluster'] = [ int(x) + 1 for x in de_tbl['cluster'] ]    

        if mg in clustering_to_k:
            k = clustering_to_k[mg]
            filename = f"{bundle_dir}/markers_{k}.tsv"
            de_tbl.to_csv(filename, sep='\t', header=True, index=False)
            manifest = _set_manifest_value(manifest, 'cluster_markers', filename, k)
        else:
            filename =  f"{bundle_dir}/{marker}.tsv"      
            de_tbl.to_csv(filename, sep='\t', header=True, index=False)
            manifest = _set_manifest_value(manifest, 'meta_markers', filename, mg)

    return manifest

def _read_file_manifest(bundle_dir):
    manifest_file=f"{bundle_dir}/MANIFEST"
    manifest=OrderedDict()

    if os.path.isfile(manifest_file):
        with open(manifest_file) as fp:
            header=fp.readline()
            for line in fp:
                line_parts=line.rstrip().split("\t")
                if len(line_parts) < 3:
                    line_parts.append('')
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

def _write_matrix_from_adata(manifest, adata, slot, bundle_dir, subdir, gene_name_field='gene_name'):
    
    layer = None
    use_raw = False

    if 'raw.' in slot:
        use_raw=True
        slot = slot.replace('raw.', '')
    
    if slot != 'X':
        layer = slot

    pathlib.Path(bundle_dir, subdir).mkdir(parents=True, exist_ok=True)    
    #ss.cmd_utils.write_mtx(adata, fname_prefix = '%s/%s/' % (bundle_dir, filename), use_raw=use_raw, use_layer=layer)
    write_mtx(adata, fname_prefix = f"{bundle_dir}/{subdir}/", use_raw=use_raw, use_layer=layer, var=[gene_name_field])
    
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
