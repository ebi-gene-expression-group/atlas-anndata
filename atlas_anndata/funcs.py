import pkg_resources
import yaml
import json
import sys
import re
import scanpy as sc
import scanpy_scripts as ss
import pandas as pd
import pathlib
import shutil
from pathlib import Path
import collections
from collections import OrderedDict
import os
import gzip
import math
import csv
import numpy as np
import glob
import scipy as sp
import sklearn as skl

from .anndata_config import (
    describe_matrices,
    describe_cellmeta,
    describe_dimreds,
    describe_genemeta,
    describe_analysis,
    load_doc,
    validate_config,
)

from .util import (
    check_slot,
    clusterings_to_ks,
    check_bundle_init,
    remove_empty_dirs,
)

from .anndata_ops import (
    update_anndata,
    derive_metadata,
    get_markers_table,
    overwrite_obs_with_magetab,
    calculate_markers,
)

from .strings import (
    schema_file,
    example_bundle_dir,
    example_manifest,
    example_config_file,
    scxa_h5ad_test,
    MISSING,
)
from jsonschema import validate


def validate_anndata_with_config(
    exp_name, bundle_dir=os.getcwd(), allow_incomplete=False
):

    """Validate an anndata against a config

    >>> config, adata = validate_anndata_with_config(
    ... 'E-MTAB-6077',
    ... example_bundle_dir
    ... ) # doctest:+ELLIPSIS +NORMALIZE_WHITESPACE
    Validating ... against .../atlas_anndata/config_schema.yaml
    Config YAML file successfully validated
    Now checking config against anndata file
    ..Checking for matrices raw.X
    ..Checking for matrices filtered
    ..Checking for matrices normalised
    ..Checking for cell_meta organism_part
    ..Checking for cell_meta louvain_resolution_0.7
    ..Checking for cell_meta louvain_resolution_1.0
    ..Checking for dimension_reductions X_umap_neighbors_n_neighbors_3
    ..Checking for dimension_reductions X_umap_neighbors_n_neighbors_10
    ..Checking for dimension_reductions X_umap_neighbors_n_neighbors_10
    ..Checking for matrices normalised
    ..Checking for gene_meta index
    ..Checking for gene_meta gene_name
    annData file successfully validated against config ...
    """

    bundle_subdir = f"{bundle_dir}/{exp_name}"
    anndata_config = f"{bundle_subdir}/anndata-config.yaml"
    config = load_doc(anndata_config)
    anndata_file = f"{bundle_subdir}/{exp_name}.project.h5ad"

    # First validate the anndata descripton file against the YAML schema

    config_status = validate_config(config, allow_incomplete=allow_incomplete)

    if config_status:
        print("Config YAML file successfully validated")
    else:
        errmsg = "Failed to validate config YAML"
        raise Exception(errmsg)

    # Now check that the things the YAML said about the annData file are true

    print("Now checking config against anndata file")
    import scanpy as sc

    adata = sc.read(anndata_file)

    # Check that data is present at the locations indicated

    for slot_type in ["matrices", "cell_meta", "dimension_reductions"]:
        if slot_type in config:
            for slot_def in config[slot_type]["entries"]:
                check_slot(
                    adata,
                    slot_type,
                    slot_def["slot"],
                    allow_incomplete=allow_incomplete,
                )

    # Check that scxa load matrix (if specified) is present

    if "load_to_scxa_db" in config["matrices"]:
        check_slot(
            adata,
            "matrices",
            config["matrices"]["load_to_scxa_db"],
            allow_incomplete=allow_incomplete,
        )

    check_slot(
        adata,
        "gene_meta",
        config["gene_meta"]["id_field"],
        allow_incomplete=allow_incomplete,
    )
    check_slot(
        adata,
        "gene_meta",
        config["gene_meta"]["name_field"],
        allow_incomplete=allow_incomplete,
    )

    if "sample_field" in config["cell_meta"]:
        check_slot(
            adata,
            "cell_meta",
            config["cell_meta"]["sample_field"],
            allow_incomplete=allow_incomplete,
        )

    # Check that some necessary version info is present

    analysis_names = [
        x["analysis"].lower() for x in config["analysis_versions"]
    ]
    required_analyses = ["reference", "mapping"]
    if not set(required_analyses).issubset(analysis_names):
        errmsg = (
            f"At least {required_analyses} must be described in"
            " analysis_versions in config file. You only have"
            f" {analysis_names}."
        )
        raise Exception(errmsg)

    print(
        f"annData file successfully validated against config {anndata_config}"
    )
    return (config, adata)


def initialise_bundle(
    exp_name,
    anndata_file,
    bundle_dir=os.getcwd(),
    atlas_style=False,
    droplet=False,
    gene_id_field="index",
    gene_name_field="gene_name",
    sample_field="sample",
    default_clustering=None,
    analysis_versions_file=None,
    **kwargs,
):

    """
    Make a yaml-format configuration file as a starting point for manual
    editing, from the content of a provided annData file.

    >>> initialise_bundle('E-MTAB-6077', scxa_h5ad_test)
    ..Checking for gene_meta gene_name
    ..Checking for gene_meta index
    """

    if anndata_file is None:
        errmsg = "No annData file supplied for bundle initialisation."
        raise Exception(errmsg)

    # Clear and create the output location

    bundle_subdir = f"{bundle_dir}/{exp_name}"

    if Path(bundle_subdir).is_dir():
        shutil.rmtree(bundle_subdir)

    pathlib.Path(bundle_subdir).mkdir(parents=True)

    # Initialise the manifest

    manifest = read_file_manifest(bundle_dir)

    # Read and describe the anndata file

    adata = sc.read(anndata_file)

    config = {
        "droplet": droplet,
        "matrices": describe_matrices(
            adata, atlas_style=atlas_style, droplet=droplet
        ),
        "cell_meta": describe_cellmeta(
            adata,
            atlas_style=atlas_style,
            droplet=droplet,
            default_clustering=default_clustering,
            sample_field=sample_field,
        ),
        "dimension_reductions": describe_dimreds(
            adata, atlas_style=atlas_style, droplet=droplet
        ),
        "gene_meta": describe_genemeta(
            adata,
            atlas_style=atlas_style,
            droplet=droplet,
            gene_id_field=gene_id_field,
            gene_name_field=gene_name_field,
        ),
        "analysis_versions": describe_analysis(
            adata,
            atlas_style=atlas_style,
            droplet=droplet,
            analysis_versions_file=analysis_versions_file,
        ),
    }

    # Write copy of annData to the bundle

    adata_filename = f"{exp_name}.project.h5ad"
    adata.write(f"{bundle_subdir}/{adata_filename}")
    set_manifest_value(manifest, "project_file", adata_filename)

    # Write the manifest

    write_file_manifest(bundle_subdir, manifest)

    # Write the config

    write_config(config, bundle_dir, exp_name)


def make_bundle_from_anndata(
    exp_name,
    step,
    bundle_dir=os.getcwd(),
    max_rank_for_stats=5,
    conda_prefix=None,
    scxa_metadata_branch="master",
    sanitize_columns=True,
    matrix_for_markers=None,
    scxa_db_scale=1000000,
    atlas_style=False,
    default_clustering=None,
    write_matrices=False,
    **kwargs,
):

    """Make the bundle from an annData file and an associated config

    >>> exp_name='E-MTAB-6077'
    >>> shutil.rmtree(exp_name, ignore_errors=True)
    >>> shutil.copytree(f"{example_bundle_dir}/{exp_name}", exp_name)
    'E-MTAB-6077'
    >>> make_bundle_from_anndata(exp_name = exp_name, step = 'final', write_matrices = True) # doctest:+ELLIPSIS +NORMALIZE_WHITESPACE
    Validating config against .../config_schema.yaml
    Config YAML file successfully validated
    Now checking config against anndata file
    ..Checking for matrices raw.X
    ..Checking for matrices filtered
    ..Checking for matrices normalised
    ..Checking for cell_meta organism_part
    ..Checking for cell_meta louvain_resolution_0.7
    ..Checking for cell_meta louvain_resolution_1.0
    ..Checking for dimension_reductions X_umap_neighbors_n_neighbors_3
    ..Checking for dimension_reductions X_umap_neighbors_n_neighbors_10
    ..Checking for dimension_reductions X_umap_neighbors_n_neighbors_10
    ..Checking for matrices normalised
    ..Checking for gene_meta index
    ..Checking for gene_meta gene_name
    annData file successfully validated against config ...
    All marker sets detailed in config present and correct
    Clearing any exisiting matrices
    Writing matrices
    .. Writing matrix from slot raw.X to subdir raw
    .. Writing matrix from slot filtered to subdir filtered
    .. normalised is a matrix we'll load into SCXA, so we'll scale it to the correct factor before we write it
    .. Writing matrix from slot normalised to subdir filtered_normalised
    Writing var(gene) metadata
    Writing obs metadata of kind: curation
    Writing obs (unsupervised clusterings)
    Writing markers and statistics
    Calculating summary stats for normalised matrix, cell groups defined by ['louvain_resolution_0.7', 'louvain_resolution_1.0']
    ..louvain_resolution_0.7
    ..louvain_resolution_1.0
    Completed summary stats calculation
    Compiling summaries for cell groupings
    ..calculating summary for cell grouping louvain_resolution_0.7
    ..limiting stats report to top 5 differential genes
    ..calculating summary for cell grouping louvain_resolution_1.0
    ..limiting stats report to top 5 differential genes
    Done compiling summaries for cell groupings
    Writing dimension reductions
    .. Writing dimension reduction from slot: X_umap_neighbors_n_neighbors_3
    .. Writing dimension reduction from slot: X_umap_neighbors_n_neighbors_10
    .. Writing dimension reduction from slot: X_umap_neighbors_n_neighbors_10
    Writing analysis versions table
    Writing annData file
    >>> shutil.rmtree(exp_name, ignore_errors=True)"""

    bundle_subdir = f"{bundle_dir}/{exp_name}"

    if step == "init":
        initialise_bundle(
            exp_name,
            atlas_style=atlas_style,
            default_clustering=default_clustering,
            **kwargs,
        )
    else:
        # Check bundle was initialied correctly
        check_bundle_init(exp_name, bundle_dir=bundle_dir)

    # Make sure the config matches the schema and anndata

    config, adata = validate_anndata_with_config(
        exp_name, bundle_dir, allow_incomplete=(step != "final")
    )

    # Read the manifest

    manifest = read_file_manifest(bundle_subdir)

    # Make any required updates to the annData object

    if step == "final":
        update_anndata(adata, config, matrix_for_markers=matrix_for_markers)

    # Write matrices at all stages after init

    if step != "init":
        write_cell_library_mapping(
            manifest=manifest,
            bundle_dir=bundle_subdir,
            adata=adata,
            config=config,
        )

    if step == "final" or write_matrices:
        write_matrices_from_adata(
            manifest=manifest,
            bundle_dir=bundle_subdir,
            adata=adata,
            config=config,
            scxa_db_scale=scxa_db_scale,
        )

    if step == "inject_magetab":

        # If curation has been done and MAGE-TAB metadata is available, then we'll
        # re-write the metadata of the object. If not, we output 'pre' magetab for curation

        adata = overwrite_obs_with_magetab(
            adata=adata,
            config=config,
            manifest=manifest,
            exp_name=exp_name,
            conda_prefix=conda_prefix,
            scxa_metadata_branch=scxa_metadata_branch,
            sanitize_columns=sanitize_columns,
            bundle_dir=bundle_subdir,
        )

        # Also add any new fields to the config
        updated_cellmeta = describe_cellmeta(
            adata,
            atlas_style=atlas_style,
            droplet=config["droplet"],
            default_clustering=default_clustering,
            sample_field=config["cell_meta"].get("sample_field"),
        )

        current_meta_slots = [
            x["slot"] for x in config["cell_meta"]["entries"]
        ]
        updated_meta_slots = [x["slot"] for x in updated_cellmeta["entries"]]
        new_meta_slots = list(
            set(updated_meta_slots) - set(current_meta_slots)
        )

        config["cell_meta"]["entries"] = config["cell_meta"]["entries"] + [
            x
            for x in updated_cellmeta["entries"]
            if x["slot"] in new_meta_slots
        ]

        set_manifest_value(
            manifest, "condensed_sdrf", f"{exp_name}.condensed-sdrf.tsv"
        )

    # Write gene metadata

    write_gene_metadata(
        manifest=manifest,
        adata=adata,
        bundle_dir=bundle_subdir,
    )

    # Write cell metadata (curated cell info)

    write_cell_metadata(
        manifest=manifest,
        adata=adata,
        bundle_dir=bundle_subdir,
        config=config,
        kind="curation",
        exp_name=exp_name,
        write_premagetab=(step == "init_magetab"),
    )

    if step != "init":

        # Write clusters (analytically derived cell groupings). For historical
        # reasons this is written differently to e.g. curated metadata

        write_clusters_from_adata(
            manifest=manifest,
            bundle_dir=bundle_subdir,
            adata=adata,
            config=config,
        )

        # Write any associated markers

        write_markers_from_adata(
            manifest=manifest,
            bundle_dir=bundle_subdir,
            adata=adata,
            config=config,
            write_marker_stats=True,
            max_rank_for_stats=max_rank_for_stats,
        )

        # Write any dim. reds from obsm

        write_obsms_from_adata(
            manifest=manifest,
            bundle_dir=bundle_subdir,
            adata=adata,
            config=config,
        )

        # Write software table

        if len(config["analysis_versions"]) > 0 and MISSING not in str(
            config["analysis_versions"]
        ):
            print("Writing analysis versions table")
            pd.DataFrame(config["analysis_versions"]).to_csv(
                f"{bundle_subdir}/software.tsv", sep="\t", index=False
            )
            set_manifest_value(
                manifest, "analysis_versions_file", "software.tsv"
            )

    # Write the anndata file back to the bundle

    print("Writing annData file")
    adata_filename = f"{exp_name}.project.h5ad"
    adata.write(f"{bundle_subdir}/{adata_filename}")
    set_manifest_value(manifest, "project_file", adata_filename)

    # Re-write the config

    write_config(config, bundle_dir, exp_name)
    set_manifest_value(
        manifest, "anndata_configuration_file", "anndata-config.yaml"
    )

    # Write the final file manifest

    if step == "final" and MISSING in str(config):
        print(
            "WARNING: Not writing manifest for bundle from incomplete"
            " config. For final bundle creation complete (or remove) all"
            f" {MISSING} entries from the config and re-run"
        )
    else:
        write_file_manifest(bundle_subdir, manifest)

    if step == "final":
        reconcile_manifest_bundle(bundle_subdir, manifest)


def write_config(config, bundle_dir, exp_name):

    anndata_config = f"{bundle_dir}/{exp_name}/anndata-config.yaml"

    with open(anndata_config, "w") as file:
        yaml.dump(config, file)


def write_cell_library_mapping(
    manifest,
    bundle_dir,
    adata,
    config,
    matrix="",
):
    """For droplet experiments, write a mapping between cell and library

    >>> test_bundle='atlas-anndata-test'
    >>> shutil.rmtree(test_bundle, ignore_errors=True)
    >>> pathlib.Path(test_bundle).mkdir()
    >>> config = load_doc(example_config_file)
    >>> manifest = read_file_manifest(example_bundle_dir)
    >>> adata = sc.read(scxa_h5ad_test)
    >>> # Our test data is not droplet. But we'll hack the structure to test this function
    >>> config['droplet'] = True
    >>> adata.obs['sample'] = adata.obs.index
    >>> config['cell_meta']['sample_field'] = 'sample'
    >>> write_cell_library_mapping(manifest, test_bundle, adata, config)
    Writing cell/library mapping for a droplet experiment
    >>> shutil.rmtree(test_bundle, ignore_errors=True)
    """

    if config["droplet"]:
        print("Writing cell/library mapping for a droplet experiment")

        # Store a cell/library mapping for each matrix. All matrices in an anndata
        # (even .raw.X, layers etc) have the same obs, so its a bit redundant to
        # write for possibly multiple matrices, but we do it to keep our pipelines
        # happy for now.

        sample_field = config["cell_meta"].get("sample_field", "sample")
        if sample_field not in adata.obs.columns:
            errmsg = (
                f"{sample_field} is not a valid obs variable, and this is a"
                " droplet experiment. This variable required to identify"
                " cells from different samples, please define it correctly in"
                " the config."
            )
            raise Exception(errmsg)

        cell_to_library = pd.DataFrame(
            {"cell": adata.obs_names, "library": adata.obs[sample_field]}
        )

        subdir = ""
        if matrix != "":
            subdir = f"matrices/{matrix}/"

        cell_to_library.to_csv(
            f"{bundle_dir}/{subdir}cell_to_library.txt",
            sep="\t",
            header=False,
            index=False,
        )
        manifest = set_manifest_value(
            manifest,
            "cell_to_library",
            f"{subdir}cell_to_library.txt",
            matrix,
        )


def write_matrices_from_adata(
    manifest,
    bundle_dir,
    adata,
    config,
    scxa_db_scale=1000000,
):
    """Given matrix slot definitions from a config file, write matrices in matrix market format

    >>> test_bundle='atlas-anndata-test'
    >>> shutil.rmtree(test_bundle, ignore_errors=True)
    >>> pathlib.Path(test_bundle).mkdir()
    >>> egconfig = load_doc(example_config_file)
    >>> adata = sc.read(scxa_h5ad_test)
    >>> write_matrices_from_adata(
    ... manifest=dict(),
    ... bundle_dir=test_bundle,
    ... adata=adata,
    ... config=egconfig)
    Clearing any exisiting matrices
    Writing matrices
    .. Writing matrix from slot raw.X to subdir raw
    .. Writing matrix from slot filtered to subdir filtered
    .. normalised is a matrix we'll load into SCXA, so we'll scale it to the correct factor before we write it
    .. Writing matrix from slot normalised to subdir normalised
    >>> shutil.rmtree(test_bundle, ignore_errors=True)
    """

    # Config tweaks can impact on matrices, so clear them all before writing

    print("Clearing any exisiting matrices")
    matrix_subdir = f"{bundle_dir}/matrices"
    if Path(matrix_subdir).is_dir():
        shutil.rmtree(matrix_subdir)

    pathlib.Path(matrix_subdir).mkdir(parents=True)

    print("Writing matrices")
    for slot_def in config["matrices"]["entries"]:

        # Skip any incomplete matrix entries.

        if MISSING in str(slot_def):
            print(f"Skipping incomplete matrix entry for {slot_def['slot']}")
            continue

        # If there's a matrix destined for the SCXA DB (which would need to be
        # count-based), apply a scaling factor if require

        if (
            "load_to_scxa_db" in config["matrices"]
            and slot_def["slot"] == config["matrices"]["load_to_scxa_db"]
        ):
            print(
                f".. {slot_def['slot']} is a matrix we'll load into SCXA, so"
                " we'll scale it to the correct factor before we write it"
            )
            scale_matrix_in_anndata(
                adata,
                slot=slot_def["slot"],
                reverse_transform=slot_def["log_transformed"],
                scale=scxa_db_scale,
            )

        write_matrix_from_adata(
            manifest=manifest,
            adata=adata,
            config=config,
            slot=slot_def["slot"],
            bundle_dir=bundle_dir,
            subdir=slot_def["name"],
            gene_name_field=config["gene_meta"]["name_field"],
        )


def write_clusters_from_adata(manifest, bundle_dir, adata, config):

    """Given obs slot definitions of clustering type from a config file, write cluster definitions to a file.

    >>> test_bundle='atlas-anndata-test'
    >>> shutil.rmtree(test_bundle, ignore_errors=True)
    >>> pathlib.Path(test_bundle).mkdir()
    >>> egconfig = load_doc(example_config_file)
    >>> adata = sc.read(scxa_h5ad_test)
    >>> write_clusters_from_adata(
    ...    manifest=dict(),
    ...    bundle_dir=test_bundle,
    ...    adata=adata,
    ...    config=egconfig)
    Writing obs (unsupervised clusterings)
    >>> shutil.rmtree(test_bundle, ignore_errors=True)
    """

    print("Writing obs (unsupervised clusterings)")

    # Find the groups in obs that correspond to clusterings
    cluster_obs = [
        x["slot"]
        for x in config["cell_meta"]["entries"]
        if x["kind"] == "clustering"
        and adata.obs[x["slot"]].dtype.name == "category"
    ]
    default_cluster_obs = [
        x["default"]
        for x in config["cell_meta"]["entries"]
        if x["kind"] == "clustering"
        and adata.obs[x["slot"]].dtype.name == "category"
    ]

    # Map clusters to k values
    clustering_to_k = clusterings_to_ks(adata, cluster_obs)

    # Create a DataFrame for cluster outputs

    clusters = adata.obs[cluster_obs].T
    clusters["K"] = [clustering_to_k[x] for x in cluster_obs]
    clusters["sel.K"] = default_cluster_obs

    clusters.to_csv(
        f"{bundle_dir}/clusters_for_bundle.txt",
        sep="\t",
        columns=["sel.K", "K"] + list(adata.obs_names),
        index=False,
    )

    set_manifest_value(
        manifest,
        "cluster_memberships",
        "clusters_for_bundle.txt",
    )


def write_obsms_from_adata(manifest, bundle_dir, adata, config):

    """Given obsm slot definitions from a config file, dimension resductions to files.

    >>> test_bundle='atlas-anndata-test'
    >>> shutil.rmtree(test_bundle, ignore_errors=True)
    >>> pathlib.Path(test_bundle).mkdir()
    >>> egconfig = load_doc(example_config_file)
    >>> adata = sc.read(scxa_h5ad_test)
    >>> write_obsms_from_adata(
    ...    manifest=dict(),
    ...    bundle_dir=test_bundle,
    ...    adata=adata,
    ...    config=egconfig)
    Writing dimension reductions
    .. Writing dimension reduction from slot: X_umap_neighbors_n_neighbors_3
    .. Writing dimension reduction from slot: X_umap_neighbors_n_neighbors_10
    .. Writing dimension reduction from slot: X_umap_neighbors_n_neighbors_10
    >>> shutil.rmtree(test_bundle, ignore_errors=True)
    """

    print("Writing dimension reductions")
    for slot_def in config["dimension_reductions"]["entries"]:
        if MISSING in str(slot_def):
            print(f"Skipping incomplete obsm entry for {slot_def['slot']}")
        else:
            write_obsm_from_adata(
                manifest,
                adata,
                obsm_name=slot_def["slot"],
                embedding_type=slot_def["kind"],
                parameterisation=json.dumps(
                    [{key: val} for key, val in slot_def["parameters"].items()]
                )
                if len(slot_def["parameters"]) > 0
                else "",
                bundle_dir=bundle_dir,
            )


# Scale a matrix to a common factor


def scale_matrix(mat, scale=1000000, reverse_transform=False):
    """
    Given a sparse matrix and a scale, derive a scaling factor from the column medians and apply scaling such that the new median matches the target. If data are log tranformed reverse that first.

    >>> adata = sc.read(scxa_h5ad_test)
    >>> print(adata.layers['normalised'][1, 2])
    7.3389006
    >>> # Scale up to per 10M
    >>> adata.layers['normalised'] = scale_matrix(adata.layers['normalised'], 10000000)
    >>> print(adata.layers['normalised'][1, 2])
    73.38901
    """
    if reverse_transform:
        mat = np.expm1(mat)
    multiplier = scale / np.median(np.ravel(mat.sum(axis=1)))
    return mat * multiplier


def scale_matrix_in_anndata(
    adata, slot, reverse_transform=False, scale=1000000
):
    """
    Given an annData object and a slot name, scale back to the standard SCXA scale using the scale_matrix() function.

    >>> adata = sc.read(scxa_h5ad_test)
    >>> print(adata.layers['normalised'][1, 2])
    7.3389006
    >>> # Scale up to per 10M
    >>> scale_matrix_in_anndata(adata, slot = 'normalised', scale = 10000000)
    >>> print(adata.layers['normalised'][1, 2])
    73.38901
    """

    if "raw." in slot:
        adata.raw.X = scale_matrix(
            adata.raw.X, scale, reverse_transform=reverse_transform
        )
    elif slot == "X":
        adata.X = scale_matrix(
            adata.X, scale, reverse_transform=reverse_transform
        )
    else:
        adata.layers[slot] = scale_matrix(
            adata.layers[slot], scale, reverse_transform=reverse_transform
        )


def write_matrix_from_adata(
    manifest,
    config,
    adata,
    slot,
    bundle_dir,
    subdir,
    gene_name_field="gene_name",
):
    """
    >>> test_bundle='atlas-anndata-test'
    >>> shutil.rmtree(test_bundle, ignore_errors=True)
    >>> pathlib.Path(test_bundle).mkdir()
    >>> adata = sc.read(scxa_h5ad_test)
    >>> egconfig = load_doc(example_config_file)
    >>> write_matrix_from_adata(
    ...     manifest=dict(),
    ...     adata=adata,
    ...     config=egconfig,
    ...     slot='X',
    ...     bundle_dir=test_bundle,
    ...     subdir='normalised',
    ...     gene_name_field = 'gene_name')
    .. Writing matrix from slot X to subdir normalised
    >>> shutil.rmtree(test_bundle, ignore_errors=True)"""

    print(f".. Writing matrix from slot {slot} to subdir {subdir}")

    layer = None
    use_raw = False

    if "raw." in slot:
        use_raw = True
        slot = slot.replace("raw.", "")

    if slot != "X":
        layer = slot

    # Make sure we have a gene name filled for every gene

    fill_gene_names(adata, gene_name_field)

    # Create the subdir

    pathlib.Path(bundle_dir, "matrices", subdir).mkdir(
        parents=True, exist_ok=True
    )
    ss.cmd_utils.write_mtx(
        adata,
        fname_prefix=f"{bundle_dir}/matrices/{subdir}/",
        use_raw=use_raw,
        use_layer=layer,
        var=[gene_name_field],
        compression={"method": "gzip"},
    )

    write_cell_library_mapping(
        manifest=manifest,
        bundle_dir=f"{bundle_dir}",
        adata=adata,
        config=config,
        matrix=subdir,
    )

    manifest = set_manifest_value(
        manifest,
        "mtx_matrix_content",
        f"matrices/{subdir}/matrix.mtx.gz",
        subdir,
    )
    manifest = set_manifest_value(
        manifest,
        "mtx_matrix_cols",
        f"matrices/{subdir}/barcodes.tsv.gz",
        subdir,
    )
    manifest = set_manifest_value(
        manifest, "mtx_matrix_rows", f"matrices/{subdir}/genes.tsv.gz", subdir
    )


def write_gene_metadata(
    manifest,
    adata,
    bundle_dir,
):
    print("Writing var(gene) metadata")

    genemeta_filename = "reference/gene_annotation.txt"
    pathlib.Path(f"{bundle_dir}/reference").mkdir(parents=True, exist_ok=True)

    nonmeta_var_patterns = [
        "mean",
        "counts",
        "n_cells",
        "highly_variable",
        "dispersion",
    ]

    nonmeta_cols = [
        x
        for x in adata.var.columns
        if any([y in x for y in nonmeta_var_patterns])
    ]
    genemeta_cols = [x for x in adata.var.columns if x not in nonmeta_cols]

    adata.var[genemeta_cols].to_csv(
        f"{bundle_dir}/{genemeta_filename}",
        sep="\t",
        header=True,
        index=True,
        index_label="gene_id",
    )
    manifest = set_manifest_value(manifest, "gene_metadata", genemeta_filename)


# Write cell metadata, including for curation as mage-tab


def write_cell_metadata(
    manifest,
    adata,
    bundle_dir,
    config,
    kind=None,
    write_premagetab=False,
    exp_name="NONAME",
):
    """Given obs slot definitions from a config file, write cell metadata to a file.

    >>> test_bundle='atlas-anndata-test'
    >>> shutil.rmtree(test_bundle, ignore_errors=True)
    >>> pathlib.Path(test_bundle).mkdir()
    >>> egconfig = load_doc(example_config_file)
    >>> adata = sc.read(scxa_h5ad_test)
    >>> write_cell_metadata(
    ...    manifest=dict(),
    ...    adata=adata,
    ...    bundle_dir=test_bundle,
    ...    config=egconfig,
    ...    kind="curation",
    ...    write_premagetab=True)
    Writing obs metadata of kind: curation
    >>> shutil.rmtree(test_bundle, ignore_errors=True)
    """

    print(f"Writing obs metadata of kind: {kind}")
    cellmeta_filename = f"{exp_name}.cell_metadata.tsv"
    presdrf_filename = f"mage-tab/{exp_name}.presdrf.txt"
    precells_filename = f"mage-tab/{exp_name}.precells.txt"

    cell_metadata, run_metadata, cell_specific_metadata = derive_metadata(
        adata, config=config, kind=None
    )

    # Output the total cell metadata anyway

    cell_metadata.to_csv(
        f"{bundle_dir}/{cellmeta_filename}",
        sep="\t",
        header=True,
        index=True,
        index_label="id",
    )
    manifest = set_manifest_value(manifest, "cell_metadata", cellmeta_filename)

    if kind == "curation" and write_premagetab:
        pathlib.Path(f"{bundle_dir}/mage-tab").mkdir(
            parents=True, exist_ok=True
        )
        if config["droplet"]:

            run_metadata.to_csv(
                f"{bundle_dir}/{presdrf_filename}",
                sep="\t",
                header=True,
                index=True,
                index_label="id",
            )

            if len(cell_specific_metadata.columns) > 0:
                cell_specific_metadata.to_csv(
                    f"{bundle_dir}/{precells_filename}",
                    sep="\t",
                    header=True,
                    index=True,
                    index_label="Cell ID",
                )
            else:
                print(
                    "Supplied anndata contained no cell-specific metadata for"
                    " this droplet experiment"
                )

        else:

            # For plate-based data the cell metadata IS the sample metadata

            cell_metadata.to_csv(
                f"{bundle_dir}/{presdrf_filename}",
                sep="\t",
                header=True,
                index=True,
                index_label="id",
            )


# Write dimension reductions from .obsm slots


def write_obsm_from_adata(
    manifest, adata, obsm_name, embedding_type, parameterisation, bundle_dir
):
    """
    Write a single dimension reduction from the specified slot

    >>> test_bundle='atlas-anndata-test'
    >>> shutil.rmtree(test_bundle, ignore_errors=True)
    >>> pathlib.Path(test_bundle).mkdir()
    >>> adata = sc.read(scxa_h5ad_test)
    >>> write_obsm_from_adata(
    ...     manifest = dict(),
    ...     adata = adata,
    ...     obsm_name='X_umap_neighbors_n_neighbors_3',
    ...     embedding_type='umap',
    ...     parameterisation='',
    ...     bundle_dir=test_bundle)
    .. Writing dimension reduction from slot: X_umap_neighbors_n_neighbors_3
    >>> shutil.rmtree(test_bundle, ignore_errors=True)
    """

    print(f".. Writing dimension reduction from slot: {obsm_name}")
    obsm_namestring = obsm_name.replace("X_", "")
    filename = f"{obsm_namestring}.tsv"
    description = f"{embedding_type}_embeddings"

    # Call the scanpy-scripts routine which writes embeddings
    ss.obj_utils.write_embedding(
        adata, key=obsm_name, embed_fn=f"{bundle_dir}/{filename}"
    )

    # Record in manifest
    manifest = set_manifest_value(
        manifest, description, filename, parameterisation
    )


# Write markers for clusterings tagged in the config


def write_markers_from_adata(
    manifest,
    bundle_dir,
    adata,
    config,
    write_marker_stats=True,
    max_rank_for_stats=5,
):

    """
    Write marker tables with p values etc from .uns

    >>> test_bundle='atlas-anndata-test'
    >>> shutil.rmtree(test_bundle, ignore_errors=True)
    >>> pathlib.Path(test_bundle).mkdir()
    >>> egconfig = load_doc(example_config_file)
    >>> adata = sc.read(scxa_h5ad_test)
    >>> write_markers_from_adata(
    ...    manifest=dict(),
    ...    bundle_dir=test_bundle,
    ...    adata=adata,
    ...    config=egconfig,
    ...    write_marker_stats=True,
    ...    max_rank_for_stats=5)
    Writing markers and statistics
    Calculating summary stats for normalised matrix, cell groups defined by ['louvain_resolution_0.7', 'louvain_resolution_1.0']
    ..louvain_resolution_0.7
    ..louvain_resolution_1.0
    Completed summary stats calculation
    Compiling summaries for cell groupings
    ..calculating summary for cell grouping louvain_resolution_0.7
    ..limiting stats report to top 5 differential genes
    ..calculating summary for cell grouping louvain_resolution_1.0
    ..limiting stats report to top 5 differential genes
    Done compiling summaries for cell groupings
    >>> shutil.rmtree(test_bundle, ignore_errors=True)
    """

    print("Writing markers and statistics")

    marker_groupings_kinds = dict(
        [
            (x["slot"], x["kind"])
            for x in config["cell_meta"]["entries"]
            if x["markers"]
        ]
    )
    marker_groupings = list(marker_groupings_kinds.keys())
    if len(marker_groupings) == 0:
        print(
            "No cell groupings have markers specified, skipping writing of"
            " markers and stats"
        )
        return

    clustering_to_k = clusterings_to_ks(adata, marker_groupings)

    # Pre-calculate the d/e tables so they can be re-used for stats

    de_tables = dict(
        zip(
            marker_groupings,
            [get_markers_table(adata, mg) for mg in marker_groupings],
        )
    )

    for cell_grouping, cell_group_kind in marker_groupings_kinds.items():

        markers_name = cell_grouping
        de_table = de_tables[cell_grouping]
        marker_type = "meta"

        if cell_group_kind == "clustering":

            # Atlas currently stores clusterings by number of clusters

            marker_type = "cluster"
            markers_name = clustering_to_k[cell_grouping]

            # Reset cluster numbering to be from 1 if required

            if de_table["cluster"].min() == "0":
                de_table["cluster"] = [int(x) + 1 for x in de_table["cluster"]]

        # Write marker table to tsv

        de_table.to_csv(
            f"{bundle_dir}/markers_{markers_name}.tsv",
            sep="\t",
            header=True,
            index=False,
        )
        manifest = set_manifest_value(
            manifest,
            f"{marker_type}_markers",
            f"markers_{markers_name}.tsv",
            markers_name,
        )

    # Now make summary statstics if we have an appropriate matrix

    if write_marker_stats and "load_to_scxa_db" in config["matrices"]:
        if MISSING in str(config["matrices"]["load_to_scxa_db"]):
            print(
                "load_to_scxa_db not set in config, can't calculate marker"
                " stats"
            )
        else:
            matrix_for_stats = config["matrices"]["load_to_scxa_db"]
            matrix_for_stats_name = [
                x["name"]
                for x in config["matrices"]["entries"]
                if x["slot"] == matrix_for_stats
            ][0]

            # Add mean and median for cell groups to the anndata

            calculate_summary_stats(
                adata,
                marker_groupings,
                matrix=matrix_for_stats,
            )

            print("Compiling summaries for cell groupings")
            marker_summary = pd.concat(
                [
                    make_markers_summary(
                        adata,
                        config["matrices"]["load_to_scxa_db"],
                        cell_grouping,
                        de_table,
                        max_rank=max_rank_for_stats,
                        cell_group_kind=marker_groupings_kinds[cell_grouping],
                    )
                    for cell_grouping, de_table in de_tables.items()
                ]
            )
            print("Done compiling summaries for cell groupings")
            statsfile = f"{matrix_for_stats_name}_stats.csv"

            marker_summary.to_csv(
                f"{bundle_dir}/{statsfile}",
                index=False,
                quoting=csv.QUOTE_NONNUMERIC,
            )
            manifest = set_manifest_value(
                manifest, "marker_stats", statsfile, matrix_for_stats_name
            )


def make_markers_summary(
    adata, layer, marker_grouping, de_table, max_rank=5, cell_group_kind=None
):
    """
    Make a table of summary statistics for a particular cell gropuing, including group-wise means.


    >>> egconfig = load_doc(example_config_file)
    >>> adata = sc.read(scxa_h5ad_test)
    >>> marker_grouping = 'louvain_resolution_0.7'
    >>> de_table = get_markers_table(adata, marker_grouping)
    >>> calculate_summary_stats(
    ...    adata,
    ...    [ marker_grouping ],
    ...    matrix='normalised')
    Calculating summary stats for normalised matrix, cell groups defined by ['louvain_resolution_0.7']
    ..louvain_resolution_0.7
    Completed summary stats calculation
    >>> make_markers_summary(
    ...    adata,
    ...    egconfig["matrices"]["load_to_scxa_db"],
    ...    marker_grouping,
    ...    de_table,
    ...    max_rank=5,
    ...    cell_group_kind='clustering')
    ..calculating summary for cell grouping louvain_resolution_0.7
    ..limiting stats report to top 5 differential genes
                    gene_id  ... median_expression
    0    ENSDARG00000017821  ...        222.825500
    0    ENSDARG00000017821  ...          0.438023
    1    ENSDARG00000079374  ...        537.408813
    1    ENSDARG00000079374  ...         35.092716
    2    ENSDARG00000045930  ...         96.161972
    2    ENSDARG00000045930  ...          0.421609
    3    ENSDARG00000076895  ...        205.787033
    3    ENSDARG00000076895  ...          2.973959
    4    ENSDARG00000070494  ...         90.376778
    4    ENSDARG00000070494  ...          0.000000
    100  ENSDARG00000030449  ...          2.710159
    100  ENSDARG00000030449  ...        556.067261
    101  ENSDARG00000051890  ...        403.266357
    101  ENSDARG00000051890  ...       2231.194580
    102  ENSDARG00000016153  ...          3.065661
    102  ENSDARG00000016153  ...        281.213318
    103  ENSDARG00000027088  ...          0.000000
    103  ENSDARG00000027088  ...        109.589203
    104  ENSDARG00000058538  ...        206.326599
    104  ENSDARG00000058538  ...       2681.499023
    <BLANKLINE>
    [20 rows x 7 columns]
    """

    print(f"..calculating summary for cell grouping {marker_grouping}")
    if max_rank:
        print(f"..limiting stats report to top {max_rank} differential genes")
        de_table = de_table[de_table["rank"] <= (max_rank - 1)]

    markers_summary = (
        de_table.drop(["ref", "scores", "logfoldchanges", "pvals"], axis=1)
        .merge(
            pd.concat(
                [
                    adata.varm[f"mean_{layer}_{marker_grouping}"]
                    .loc[de_table["genes"].unique()]
                    .melt(ignore_index=False),
                    adata.varm[f"median_{layer}_{marker_grouping}"]
                    .loc[de_table["genes"].unique()]
                    .melt(ignore_index=False),
                ],
                axis=1,
            )
            .iloc[:, [0, 1, 3]]
            .set_axis(
                ["cluster_id", "mean_expression", "median_expression"], axis=1
            ),
            left_on="genes",
            right_index=True,
        )
        .rename(
            columns={
                "genes": "gene_id",
                "cluster": "group_where_marker",
                "pvals_adj": "marker_p_value",
            }
        )
    )

    # For unsupervised clusterings, record the grouping as k and increment the
    # group numbers so they start from 1 if they're integers

    if cell_group_kind == "clustering":
        markers_summary["grouping_where_marker"] = len(
            adata.obs[marker_grouping].unique()
        )
        if (
            markers_summary["cluster_id"].str.isnumeric().all()
            and min([int(x) for x in markers_summary["cluster_id"]]) == 0
        ):
            markers_summary["cluster_id"] = [
                int(x) + 1 for x in markers_summary["cluster_id"]
            ]
    else:
        markers_summary["grouping_where_marker"] = marker_grouping

    # Convert cell group columns to strings

    for col in ["grouping_where_marker", "group_where_marker", "cluster_id"]:
        markers_summary[col] = markers_summary[col].astype(str)

    return markers_summary[
        [
            "gene_id",
            "grouping_where_marker",
            "group_where_marker",
            "cluster_id",
            "marker_p_value",
            "mean_expression",
            "median_expression",
        ]
    ]


def read_file_manifest(bundle_dir=None, manifest_file=None):
    """
    Read a manifest into a dictionary, or initialise an empty one

    >>> read_file_manifest('test_bundle', manifest_file = example_manifest) # doctest:+ELLIPSIS +NORMALIZE_WHITESPACE
    OrderedDict([('cell_metadata', OrderedDict([('', 'E-MTAB-6077.cell_metadata.tsv')])), ...
    """

    manifest_file = (
        f"{bundle_dir}/MANIFEST" if manifest_file is None else manifest_file
    )
    manifest = OrderedDict()

    if os.path.isfile(manifest_file):
        with open(manifest_file) as fp:
            _ = fp.readline()
            for line in fp:
                line_parts = line.rstrip().split("\t")
                if len(line_parts) < 3:
                    line_parts.append("")
                manifest = set_manifest_value(
                    manifest, line_parts[0], line_parts[1], line_parts[2]
                )

    return manifest


def reconcile_manifest_bundle(bundle_dir, manifest):

    """
    Remove files in a bundle directory which are absent from the manifest

    >>> exp_name='E-MTAB-6077'
    >>> shutil.rmtree(exp_name, ignore_errors=True)
    >>> shutil.copytree(f"{example_bundle_dir}/{exp_name}", exp_name)
    'E-MTAB-6077'
    >>> manifest = read_file_manifest(manifest_file = f"{exp_name}/MANIFEST")
    >>> reconcile_manifest_bundle(exp_name, manifest)
    >>> shutil.rmtree(exp_name, ignore_errors=True)
    """

    cwd = os.getcwd()
    manifest_files = [
        item
        for sublist in [list(x.values()) for x in manifest.values()]
        for item in sublist
    ]
    os.chdir(bundle_dir)
    bundle_files = [
        x for x in glob.glob("**", recursive=True) if not os.path.isdir(x)
    ]

    # Remove any files not in the manifest

    for non_manifest_file in [
        x for x in bundle_files if x not in manifest_files
    ]:
        if non_manifest_file != "MANIFEST":
            print(f"Cleaning up {non_manifest_file}")
            os.remove(non_manifest_file)
    os.chdir(cwd)

    # Finally clean up any empty directories

    remove_empty_dirs(bundle_dir)


def write_file_manifest(bundle_dir, manifest):
    """
    Write a manifest into a dictionary

    >>> test_bundle='atlas-anndata-test'
    >>> shutil.rmtree(test_bundle, ignore_errors=True)
    >>> pathlib.Path(test_bundle).mkdir()
    >>> manifest = read_file_manifest(manifest_file = example_manifest)
    >>> write_file_manifest(test_bundle, manifest)
    >>> shutil.rmtree(test_bundle, ignore_errors=True)
    """

    manifest_file = f"{bundle_dir}/MANIFEST"

    with open(manifest_file, "w") as fh:
        fh.write("Description\tFile\tParameterisation\n")

        for description, v in manifest.items():
            for parameterisation, filename in v.items():
                fh.write(f"{description}\t{filename}\t{parameterisation}\n")


def set_manifest_value(manifest, description, filename, parameterisation=""):
    """
    Add a manifest value

    >>> manifest = read_file_manifest(example_manifest)
    >>> set_manifest_value(
    ...    manifest, "analysis_versions_file", f"foo.tsv")
    OrderedDict([('analysis_versions_file', OrderedDict([('', 'foo.tsv')]))])
    """

    if description not in manifest:
        manifest[description] = OrderedDict()

    manifest[description][parameterisation] = filename

    return manifest


def fill_gene_names(adata, gene_name_field="gene_name"):
    """Make sure that a gene name field is filled for all genes by filling blanks with gene IDs

    >>> gene_name_field = 'gene_name'
    >>> adata = sc.read(scxa_h5ad_test)
    >>> genes_with_missing_names = len(list(
    ...    adata.var_names[pd.isnull(adata.var[gene_name_field])]))
    >>> print(f"Starting with {genes_with_missing_names} genes without names")
    Starting with 48 genes without names
    >>> adata = fill_gene_names(adata, gene_name_field=gene_name_field)
    >>> genes_with_missing_names = len(list(
    ...    adata.var_names[pd.isnull(adata.var[gene_name_field])]))
    >>> print(f"Ending with {genes_with_missing_names} genes without names")
    Ending with 0 genes without names
    """

    if gene_name_field in adata.var.columns:

        # For any genes without names, assign the ID
        genes_with_missing_names = list(
            adata.var_names[pd.isnull(adata.var[gene_name_field])]
        )
        adata.var[gene_name_field] = adata.var[
            gene_name_field
        ].cat.add_categories(genes_with_missing_names)
        adata.var.loc[
            genes_with_missing_names, gene_name_field
        ] = genes_with_missing_names
    else:
        adata.var[gene_name_field] = adata.var_names

    return adata


def calculate_summary_stats(adata, obs, matrix="normalised"):

    """Calculate mean and median for cell groups defined by a variable in .obs

    >>> adata = sc.read(scxa_h5ad_test)
    >>> marker_grouping = 'louvain_resolution_0.7'
    >>> de_table = get_markers_table(adata, marker_grouping)
    >>> calculate_summary_stats(
    ...    adata,
    ...    [ marker_grouping ],
    ...    matrix='normalised')
    Calculating summary stats for normalised matrix, cell groups defined by ['louvain_resolution_0.7']
    ..louvain_resolution_0.7
    Completed summary stats calculation
    >>> adata.varm['mean_normalised_louvain_resolution_0.7'] # doctest:+ELLIPSIS +NORMALIZE_WHITESPACE
                                 0           1
    ENSDARG00000000001    5.210635    4.868471
    ENSDARG00000000002    5.107658    0.254911
    ENSDARG00000000018  240.056006  319.279625
    ENSDARG00000000019   34.447536    9.438827
    ENSDARG00000000068   33.625522    1.097745
    ...
    >>> adata.varm['median_normalised_louvain_resolution_0.7'] # doctest:+ELLIPSIS +NORMALIZE_WHITESPACE
                             0           1
    ENSDARG00000000001    0.000000    0.000000
    ENSDARG00000000002    0.000000    0.000000
    ENSDARG00000000018  141.642776  159.815002
    ENSDARG00000000019    0.000000    0.000000
    ENSDARG00000000068    0.000000    0.000000
    ...
    """

    print(
        f"Calculating summary stats for {matrix} matrix, cell groups defined"
        f" by {obs}"
    )

    def matrix_summary(matrix, summary_type):
        if summary_type == "mean":
            return matrix.mean(axis=0, dtype=np.float64)
        else:
            return skl.utils.sparsefuncs.csc_median_axis_0(
                sp.sparse.csc_matrix(matrix, dtype=np.float64)
            )

    for ob in obs:
        print(f"..{ob}")
        layer = None
        use_raw = False

        if matrix == "raw.x":
            use_raw = True

        elif matrix in adata.layers:
            layer = matrix

        for summary_type in ["mean", "median"]:
            adata.varm[f"{summary_type}_{matrix}_{ob}"] = pd.DataFrame(
                (
                    np.ravel(
                        matrix_summary(
                            sc.get._get_obs_rep(
                                adata[adata.obs[ob] == group],
                                layer=layer,
                                use_raw=use_raw,
                            ),
                            summary_type,
                        )
                    )
                    for group in adata.obs[ob].cat.categories
                ),
                columns=adata.var_names,
                index=adata.obs[ob].cat.categories,
            ).transpose()

    print("Completed summary stats calculation")
