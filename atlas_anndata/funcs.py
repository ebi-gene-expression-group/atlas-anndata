import pkg_resources
import jsonschema
from jsonschema import validate
import yaml
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

schema_file = pkg_resources.resource_filename(
    "atlas_anndata", "config_schema.yaml"
)
example_config_file = pkg_resources.resource_filename(
    "atlas_anndata", "example_config.yaml"
)
scxa_h5ad_test = pkg_resources.resource_filename(
    "atlas_anndata", "data/E-MTAB-6077.project.h5ad"
)


def load_doc(filename):
    with open(filename, "r") as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)


def validate_config(config):

    """Validate a config against our schema


    >>> egconfig = load_doc(example_config_file)
    >>> validate_config(egconfig)
    True
    """

    schema = load_doc(schema_file)

    try:
        validate(instance=config, schema=schema)
    except jsonschema.exceptions.ValidationError as err:
        print(err, file=sys.stderr)
        return False
    return True


def check_slot(adata, slot_type, slot_name):

    """Check for a slot in an anndata object

    >>> adata = sc.read(scxa_h5ad_test)
    >>> check = check_slot(adata, 'matrices', 'X')
    Checking for matrices X
    >>> check
    True
    """

    print(f"Checking for {slot_type} {slot_name}")

    if slot_type == "matrices":
        if slot_name == "X":
            return hasattr(adata, "X")
        elif slot_name == "raw.X":
            return hasattr(adata.raw, "X")
        else:
            return slot_name in adata.layers

    elif slot_type == "dimension_reductions":
        return slot_name in adata.obsm

    elif slot_type == "cell_groups":
        return slot_name in adata.obs.columns

    else:
        print(f"{slot_type} slot type not recognised", file=sys.stderr)
        return False


def validate_anndata_with_config(config_file, anndata_file):

    """Validate an anndata against a config

    >>> config, adata = validate_anndata_with_config(
    ... example_config_file,
    ... scxa_h5ad_test
    ... ) # doctest:+ELLIPSIS +NORMALIZE_WHITESPACE
    Validating .../atlas_anndata/example_config.yaml against
    .../atlas_anndata/config_schema.yaml
    Config YAML file successfully validated
    Now checking config against anndata file
    Checking for matrices raw.X
    Checking for matrices filtered
    Checking for matrices normalised
    Checking for cell_groups organism_part
    Checking for cell_groups louvain_resolution_0.7
    Checking for cell_groups louvain_resolution_1.0
    Checking for dimension_reductions X_umap_neighbors_n_neighbors_3
    Checking for dimension_reductions X_umap_neighbors_n_neighbors_10
    Checking for dimension_reductions X_umap_neighbors_n_neighbors_10
    annData file successfully validated against config ...
    """

    config = load_doc(config_file)

    # First validate the anndata descripton file against the YAML schema

    print(f"Validating {config_file} against {schema_file}")
    config_status = validate_config(config)

    if config_status:
        print("Config YAML file successfully validated")
    else:
        print("FAILURE")
        sys.exit(1)

    # Now check that the things the YAML said about the annData file are true

    print("Now checking config against anndata file")
    import scanpy as sc

    adata = sc.read(anndata_file)

    # Check that data is present at the locations indicated

    for slot_type in ["matrices", "cell_groups", "dimension_reductions"]:
        if slot_type in config:
            for slot_def in config[slot_type]["entries"]:
                if not check_slot(adata, slot_type, slot_def["slot"]):
                    print(
                        f"{slot_type} entry {slot_def['slot']} not present in"
                        f" anndata file {anndata_file}"
                    )
                    sys.exit(1)

    print(f"annData file successfully validated against config {config_file}")
    return (config, adata)


def obs_markers(adata, obs):

    """
    Check if a particular cell metadata field has an associated marker set

    >>> adata = sc.read(scxa_h5ad_test)
    >>> obs_markers(adata, 'louvain_resolution_1.0')
    'markers_louvain_resolution_1.0'
    """

    markers_slot = f"markers_{obs}"

    if markers_slot in adata.uns.keys():
        return markers_slot
    else:
        return False


def slot_kind_from_name(slot_type, slot_name):

    """
    Try to infer the kind/sub-type of a slot by comparing to known patterns


    >>> slot_kind_from_name('dimension_reductions', 'X_tsne_blah')
    'tsne'
    """

    from atlas_anndata import MISSING_STRING

    kind = MISSING_STRING
    search_map = {}

    if slot_type == "dimension_reductions":
        kind = f"{kind}: 'pca', 'tsne' or 'umap'"
        search_map = {".*pca": "pca", ".*umap": "umap", ".*tsne": "tsne"}

    elif slot_type == "cell_groups":
        kind = "curation"
        search_map = {
            "leiden": "clustering",
            "cluster": "clustering",
            "log1p": "analysis",
            "^n_": "analysis",
            "total_": "analysis",
            "pct_": "analysis",
            "mito": "analysis",
            "louvain": "clustering",
        }
    elif slot_type == "gene_meta":
        kind = "curation"
        search_map = {
            "mean": "analysis",
            "counts": "analysis",
            "n_cells": "analysis",
            "highly_variable": "analysis",
            "dispersion": "analysis",
        }

    for kind_pattern in search_map.keys():
        if re.match(r"{}".format(kind_pattern), slot_name.lower()):
            kind = search_map[kind_pattern]
            break

    return kind


def extract_parameterisation(slot_type, slot_name, atlas_style=False):

    """
    For annData objects from Single Cell Expression Atlas, infer paramerisation
    from slot naming.

    >>> extract_parameterisation(
    ... 'cell_groups',
    ... 'louvain_resolution_1.0',
    ... atlas_style = True )
    {'resolution': 1.0}
    """

    parameters = {}

    if atlas_style:

        if slot_type == "dimension_reductions":
            m = re.search(
                r"X_(umap|tsne|pca)_(.*)_(.*)",
                slot_name.replace("umap_neighbors", "umap"),
            )
            if m:
                parameters[m.group(2)] = string_to_numeric(m.group(3))

        elif slot_type == "cell_groups":
            m = re.search(r"(louvain|leiden)_(.*)_(.*)", slot_name)
            if m:
                parameters[m.group(2)] = string_to_numeric(m.group(3))

    return parameters


def string_to_numeric(numberstring):

    """
    Convert strings to int or float

    >>> string_to_numeric('1.0')
    1.0
    """

    if numberstring.isdigit():
        return int(numberstring)
    else:
        return float(numberstring)


def make_starting_config_from_anndata(
    anndata_file, config_file, atlas_style=False, exp_name=None, droplet=False
):

    """
    Make a yaml-format configuration file as a starting point for manual
    editing, from the content of a provided annData file.

    >>> make_starting_config_from_anndata(scxa_h5ad_test, '/tmp/foo.yaml')
    """

    adata = sc.read(anndata_file)
    from atlas_anndata import MISSING_STRING

    config = {
        "exp-name": exp_name,
        "droplet": droplet,
        "matrices": {"load_to_scxa_db": None, "entries": []},
        "cell_groups": {"entries": []},
        "dimension_reductions": {"entries": []},
    }

    # Describe the matrices

    matrix_slots = []

    if hasattr(adata, "raw") and hasattr(adata.raw, "X"):
        matrix_slots.append("raw.X")

    matrix_slots = matrix_slots + list(adata.layers.keys())
    matrix_slots.append("X")

    if atlas_style and "normalised" in matrix_slots:
        config["matrices"]["load_to_scxa_db"] = "normalised"

    for ms in matrix_slots:
        matrix_entry = {
            "slot": ms,
            "name": ms if ms in adata.layers else MISSING_STRING,
            "measure": "counts" if atlas_style else MISSING_STRING,
            "cell_filtered": True,
            "gene_filtered": False,
            "normalised": False,
            "log_transformed": False,
            "scaled": True,
            "parameters": {},
        }

        # Use assumptions we can rely on for Atlas

        if atlas_style:
            if ms in ["filtered", "normalised", "X"]:
                matrix_entry["gene_filtered"] = True
            if ms in ["normalised", "X"]:
                matrix_entry["normalised"] = True
            if ms == "X":
                matrix_entry["log_transformed"] = True
                if droplet:
                    matrix_entry["scaled"] = True

        config["matrices"]["entries"].append(matrix_entry)

    # Describe cell-wise metadata columns

    for obs in adata.obs.columns:

        obs_entry = {
            "slot": obs,
            "kind": slot_kind_from_name("cell_groups", obs),
            "parameters": extract_parameterisation(
                "cell_groups", obs, atlas_style
            ),
            "default": False,
        }

        markers_slot = obs_markers(adata, obs)
        if markers_slot:
            obs_entry["markers"] = True
            obs_entry["markers_slot"] = markers_slot
        else:
            obs_entry["markers"] = False

        config["cell_groups"]["entries"].append(obs_entry)

    # Find the groups in obs that correspond to clusterings
    config_obs = [x["slot"] for x in config["cell_groups"]["entries"]]
    cluster_obs = [
        x["slot"]
        for x in config["cell_groups"]["entries"]
        if x["kind"] == "clustering"
    ]

    cluster_obs = list(
        select_clusterings(
            adata, clusters=cluster_obs, atlas_style=atlas_style
        ).keys()
    )

    # Try to set a default clustering
    if atlas_style and "louvain_resolution_1.0" in cluster_obs:
        default_clustering = "louvain_resolution_1.0"
    else:
        default_clustering = cluster_obs[
            math.floor((len(cluster_obs) - 1) / 2)
        ]

    config["cell_groups"]["entries"][config_obs.index(default_clustering)][
        "default"
    ] = True

    # Describe dimension reductions stored in .obsm

    for obsm in adata.obsm.keys():
        obsm_entry = {
            "slot": obsm,
            "kind": slot_kind_from_name("dimension_reductions", obsm),
            "parameters": extract_parameterisation(
                "dimension_reductions", obsm, atlas_style
            ),
        }

        config["dimension_reductions"]["entries"].append(obsm_entry)

    with open(config_file, "w") as file:
        yaml.dump(config, file)


def make_bundle_from_anndata(
    anndata_file,
    anndata_description_yaml,
    bundle_dir,
    gene_name_field=None,
    **kwargs,
):
    # Make sure the config matches the schema and anndata

    config, adata = validate_anndata_with_config(
        anndata_description_yaml, anndata_file
    )

    # Clear and create the output location

    if Path(bundle_dir).is_dir():
        shutil.rmtree(bundle_dir)

    pathlib.Path(f"{bundle_dir}").mkdir(parents=True)

    # Initialise the manifest

    manifest = read_file_manifest(bundle_dir)

    # Write matrices

    print("Writing matrices")
    write_matrices_from_adata(
        manifest=manifest,
        bundle_dir=bundle_dir,
        adata=adata,
        config=config,
        gene_name_field=gene_name_field,
    )

    # Write clusters (analytically derived cell groupings). For historical
    # reasons this is written differently to e.g. curated metadata

    print("Writing obs (unsupervised clusterings)")
    write_clusters_from_adata(
        manifest=manifest,
        bundle_dir=bundle_dir,
        adata=adata,
        config=config,
    )

    # Write cell metadata (curated cell info)

    print("Writing obs (curated metadata)")
    write_cell_metadata(
        manifest=manifest,
        adata=adata,
        bundle_dir=bundle_dir,
        config=config,
        kind="curation",
        write_premagetab=True,
    )

    # Write any associated markers

    print("Writing markers and statistics")
    write_markers_from_adata(
        manifest=manifest,
        bundle_dir=bundle_dir,
        adata=adata,
        config=config,
        write_marker_stats=True,
    )

    # Write any dim. reds from obsm

    print("Writing dimension reductions")
    write_obsms_from_adata(
        manifest=manifest,
        bundle_dir=bundle_dir,
        adata=adata,
        config=config,
    )

    print("Writing annData file")
    adata_filename = f"{config['exp-name']}.project.h5ad"
    adata.write(f"{bundle_dir}/{adata_filename}")
    set_manifest_value(manifest, "project_file", adata_filename)

    # Write the final file manifest

    write_file_manifest(bundle_dir, manifest)


def write_matrices_from_adata(
    manifest, bundle_dir, adata, config, gene_name_field
):
    for slot_def in config["matrices"]["entries"]:
        write_matrix_from_adata(
            manifest=manifest,
            adata=adata,
            slot=slot_def["slot"],
            bundle_dir=bundle_dir,
            subdir=slot_def["name"],
            gene_name_field=gene_name_field,
        )


def write_clusters_from_adata(manifest, bundle_dir, adata, config):

    # Find the groups in obs that correspond to clusterings
    cluster_obs = [
        x["slot"]
        for x in config["cell_groups"]["entries"]
        if x["kind"] == "clustering"
    ]
    default_cluster_obs = [
        x["default"]
        for x in config["cell_groups"]["entries"]
        if x["kind"] == "clustering"
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

    for slot_def in config["dimension_reductions"]["entries"]:
        write_obsm_from_adata(
            manifest,
            adata,
            obsm_name=slot_def["slot"],
            embedding_type=slot_def["kind"],
            parameterisation=list(slot_def["parameters"].values())[0]
            if len(slot_def["parameters"]) > 0
            else "",
            bundle_dir=bundle_dir,
        )


def write_matrix_from_adata(
    manifest, adata, slot, bundle_dir, subdir, gene_name_field="gene_name"
):

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

    pathlib.Path(bundle_dir, subdir).mkdir(parents=True, exist_ok=True)
    write_mtx(
        adata,
        fname_prefix=f"{bundle_dir}/{subdir}/",
        use_raw=use_raw,
        use_layer=layer,
        var=[gene_name_field],
    )

    for filename in ["matrix.mtx", "barcodes.tsv", "genes.tsv"]:
        subfile = f"{subdir}/{filename}"
        filepath = f"{bundle_dir}/{subfile}"
        with open(filepath, "rb") as f_in, gzip.open(
            "%s.gz" % filepath, "wb"
        ) as f_out:
            f_out.writelines(f_in)
        os.remove(filepath)

    manifest = set_manifest_value(
        manifest, "mtx_matrix_content", f"{subdir}/matrix.mtx.gz", subdir
    )
    manifest = set_manifest_value(
        manifest, "mtx_matrix_cols", f"{subdir}/barcodes.tsv.gz", subdir
    )
    manifest = set_manifest_value(
        manifest, "mtx_matrix_rows", f"{subdir}/genes.mtx.gz", subdir
    )

    return manifest


def create_bundle(
    anndata_file,
    anndata_description_yaml,
    bundle_dir,
    write_cellmeta=True,
    write_genemeta=True,
    nonmeta_var_patterns=None,
    nonmeta_obs_patterns=None,
    exp_desc=None,
    software_versions_file=None,
    raw_matrix_slot=None,
    filtered_matrix_slot=None,
    normalised_matrix_slot=None,
    final_matrix_slot=None,
    write_obsms=True,
    obsms=None,
    write_clusters=True,
    clusters=None,
    clusters_field_pattern="louvain",
    default_clustering=None,
    write_markers=True,
    marker_clusterings=None,
    metadata_marker_fields=None,
    write_marker_stats=True,
    marker_stats_layers=None,
    max_rank_for_stats=4,
    atlas_style=True,
    write_anndata=True,
):
    """Build a bundle directory compatible with Single Cell Expression Atlas
    (SCXA) build proceseses

    \b
    anndata_file - A file of the annData hdf5 specification, with all necessary
                   information for SCXA.
    bundle_dir   - A directory in which to create the bundle.
    """

    adata = sc.read(anndata_file)
    config = load_doc(anndata_description_yaml)

    # Remove anything completely empty in .obs or .var

    # Add anything we need to augment the adata

    if software_versions_file is not None:
        adata.uns["software_versions"] = pd.read_csv(
            software_versions_file, sep="\t"
        )

    if Path(bundle_dir).is_dir():
        shutil.rmtree(bundle_dir)

    pathlib.Path(f"{bundle_dir}/mage-tab").mkdir(parents=True)
    pathlib.Path(f"{bundle_dir}/reference").mkdir()

    manifest = read_file_manifest(bundle_dir)

    if write_cellmeta:
        print("Writing cell and run/sample metadata")
        manifest = write_cell_metadata(
            manifest,
            adata,
            bundle_dir,
            config=config,
            nonmeta_obs_patterns=nonmeta_obs_patterns,
        )

    if write_genemeta:
        print("Writing gene-wise metadata")

        genemeta_filename = f"{bundle_dir}/reference/gene_annotation.txt"
        if nonmeta_var_patterns is None:
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
            genemeta_filename,
            sep="\t",
            header=True,
            index=True,
            index_label="gene_id",
        )
        manifest = set_manifest_value(
            manifest, "gene_metadata", "reference/gene_annotation.txt"
        )

    if write_obsms:
        print("Writing obsms to file")
        obsms = adata.obsm.keys() if obsms is None else obsms

        for obsm in obsms:
            print(f"Writing {obsm}")
            manifest = write_obsm_from_adata(
                manifest,
                adata,
                obsm_name=obsm,
                bundle_dir=bundle_dir,
                atlas_style=atlas_style,
            )

    if write_markers:
        print("Writing markers to file")
        manifest = write_markers_from_adata(
            manifest,
            adata,
            clusters=clusters,
            marker_clusterings=marker_clusterings,
            metadata_marker_fields=metadata_marker_fields,
            bundle_dir=bundle_dir,
            atlas_style=atlas_style,
            max_rank_for_stats=max_rank_for_stats,
            marker_stats_layers=marker_stats_layers,
            write_marker_stats=write_marker_stats,
        )

    if "software_versions" in adata.uns:
        print("Writing provided software info to bundle")
        software_versions_outfile = f"{bundle_dir}/software.tsv"
        adata.uns["software_versions"].to_csv(
            software_versions_outfile, sep="\t", header=True, index=False
        )
        manifest = set_manifest_value(
            manifest, "software_versions_file", "software.tsv"
        )

    # Write anndata

    if write_anndata:

        print("Writing annData file to bundle")
        adata_filename = f"{config['exp-name']}.project.h5ad"
        adata.write(f"{bundle_dir}/{adata_filename}")
        manifest = set_manifest_value(manifest, "project_file", adata_filename)

    # Write the final file manifest

    write_file_manifest(bundle_dir, manifest)


# Write cell metadata, including for curation as mage-tab


def write_cell_metadata(
    manifest, adata, bundle_dir, config, kind=None, write_premagetab=False
):

    # By default print all obs columns, but that's probably not we want in most
    # cases because of mixture of data types there (from curation, QC,
    # clustering etc)

    if kind is None:
        obs_columns = list(adata.obs.columns)
    else:
        obs_columns = [
            slot_def["slot"]
            for slot_def in config["cell_groups"]["entries"]
            if slot_def["kind"] == kind
        ]

    print("Writing cell metdata to be used in curation")
    cellmeta_filename = f"{config['exp-name']}.cell_metadata.tsv"
    presdrf_filename = f"mage-tab/{config['exp-name']}.presdrf.txt"
    precells_filename = f"mage-tab/{config['exp-name']}.precells.txt"
    pathlib.Path(f"{bundle_dir}/mage-tab").mkdir(parents=True, exist_ok=True)

    cell_metadata = adata.obs[obs_columns].copy()
    cell_metadata.to_csv(
        f"{bundle_dir}/{cellmeta_filename}",
        sep="\t",
        header=True,
        index=True,
        index_label="id",
    )
    manifest = set_manifest_value(manifest, "cell_metadata", cellmeta_filename)

    if kind == "curation" and write_premagetab:
        if config["droplet"]:

            # Split cell IDs to runs and barcodes
            try:
                runs, barcodes = zip(*(s.split("-") for s in adata.obs_names))
            except ValueError as e:
                print(
                    f"Error deriving run and barcode lists: {e}. Cell names"
                    " likely don't match expected <run or sample>_<barcode>"
                    " naming format for droplet experiments. First cell name"
                    f" is: {adata.obs_names[0]}"
                )
                sys.exit(1)

            cell_metadata["cell barcode"] = barcodes

            # Split cell metadata by run ID and create run-wise metadata with
            # any invariant value across all cells of a run

            unique_runs = list(set(runs))
            submetas = [
                cell_metadata[[y == x for y in runs]] for x in unique_runs
            ]
            run_meta = pd.concat(
                [
                    df[[x for x in df.columns if len(set(df[x])) == 1]].head(1)
                    for df in submetas
                ],
                join="inner",
            )
            run_meta["run"] = unique_runs
            run_meta.set_index("run", inplace=True)

            run_meta.to_csv(
                f"{bundle_dir}/{presdrf_filename}",
                sep="\t",
                header=True,
                index=True,
                index_label="id",
            )
            cell_specific_metadata = cell_metadata[
                [
                    x
                    for x in cell_metadata.columns
                    if x not in obs_columns + list(run_meta.columns)
                ]
            ]

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


def clusterings_to_ks(adata, obs_names):
    return dict(
        zip(obs_names, [len(adata.obs[c].unique()) for c in obs_names])
    )


def select_clusterings(adata, clusters, atlas_style=True):

    clusterings = list(
        dict(
            sorted(
                dict(
                    zip(
                        clusters,
                        [abs(1 - float(c.split("_")[-1])) for c in clusters],
                    )
                ).items(),
                key=lambda x: x[1],
            )
        ).keys()
    )

    clustering_to_nclust = clusterings_to_ks(adata, clusterings)

    if atlas_style:

        # Only keep the first marker set for a given k. For Atlas, ranked as
        # above, this will be the makers from the resolution closest to 1 of
        # clashing sets.

        kept_clusterings = {}
        for k, v in clustering_to_nclust.items():
            if v not in kept_clusterings.values():
                kept_clusterings[k] = v

        clustering_to_nclust = kept_clusterings

    return dict(sorted(clustering_to_nclust.items(), key=lambda item: item[1]))


# Write markers for clusterings tagged in the config


def write_markers_from_adata(
    manifest,
    bundle_dir,
    adata,
    config,
    write_marker_stats=True,
    max_rank_for_stats=4,
):
    marker_groupings_kinds = [
        (x["slot"], x["kind"])
        for x in config["cell_groups"]["entries"]
        if x["markers"]
    ]
    marker_groupings = [x[0] for x in marker_groupings_kinds]
    clustering_to_k = clusterings_to_ks(adata, marker_groupings)

    # Pre-calculate the d/e tables so they can be re-used for stats

    de_tables = dict(
        zip(
            marker_groupings,
            [get_markers_table(adata, mg) for mg in marker_groupings],
        )
    )

    for cell_grouping, cell_group_kind in marker_groupings_kinds:

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

        matrix_for_stats = config["matrices"]["load_to_scxa_db"]

        # Add mean and median for cell groups to the anndata

        calculate_summary_stats(
            adata,
            marker_groupings,
            matrix=matrix_for_stats,
        )

        marker_summary = pd.concat(
            [
                make_markers_summary(
                    adata,
                    config["matrices"]["load_to_scxa_db"],
                    cell_grouping,
                    de_table,
                    max_rank=max_rank_for_stats,
                    cell_group_kind=cell_group_kind,
                )
                for mg, de_table in de_tables.items()
            ]
        )
        statsfile = f"{bundle_dir}/{matrix_for_stats}_stats.csv"

        marker_summary.to_csv(statsfile, index=False)
        manifest = set_manifest_value(
            manifest, "marker_stats", statsfile, matrix_for_stats
        )


def get_markers_table(adata, marker_grouping):

    de_table = ss.lib._diffexp.extract_de_table(
        adata.uns[f"markers_{marker_grouping}"]
    )
    de_table = de_table.loc[de_table.genes.astype(str) != "nan", :]

    return de_table


def make_markers_summary(
    adata, layer, marker_grouping, de_table, max_rank=4, cell_group_kind=None
):

    summary_stats = (
        pd.concat(
            [
                adata.varm[f"mean_{layer}_{marker_grouping}"].melt(
                    ignore_index=False
                ),
                adata.varm[f"median_{layer}_{marker_grouping}"].melt(
                    ignore_index=False
                ),
            ],
            axis=1,
        )
        .iloc[:, [0, 1, 3]]
        .set_axis(
            ["cluster_id", "mean_expression", "median_expression"], axis=1
        )
    )

    new_colnames = {
        "genes": "gene_id",
        "cluster": "group_where_marker",
        "pvals_adj": "marker_p_value",
    }

    markers_summary = (
        de_table.merge(summary_stats, left_on="genes", right_index=True)
        .drop(["ref", "scores", "logfoldchanges", "pvals"], axis=1)
        .rename(columns=new_colnames)
    )

    if max_rank:
        markers_summary = markers_summary[markers_summary["rank"] <= max_rank]

    # For unsupervised clusterings, record the grouping as k and increment the
    # group numbers so they start from 1

    if cell_group_kind == "clustering":
        markers_summary["grouping_where_marker"] = len(
            adata.obs[marker_grouping].unique()
        )
        if min([int(x) for x in markers_summary["cluster_id"]]) == 0:
            markers_summary["cluster_id"] = [
                int(x) + 1 for x in markers_summary["cluster_id"]
            ]

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


def read_file_manifest(bundle_dir):
    manifest_file = f"{bundle_dir}/MANIFEST"
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


def write_file_manifest(bundle_dir, manifest):
    manifest_file = f"{bundle_dir}/MANIFEST"

    with open(manifest_file, "w") as fh:
        fh.write("Description\tFile\tParameterisation\n")

        for description, v in manifest.items():
            for parameterisation, filename in v.items():
                fh.write(f"{description}\t{filename}\t{parameterisation}\n")


def set_manifest_value(manifest, description, filename, parameterisation=""):

    if description not in manifest:
        manifest[description] = OrderedDict()

    manifest[description][parameterisation] = filename

    return manifest


def fill_gene_names(adata, gene_name_field="gene_name"):

    if gene_name_field not in adata.var_names:

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


def write_mtx(
    adata, fname_prefix="", var=None, obs=None, use_raw=False, use_layer=None
):
    """Export AnnData object to mtx formt
    * Parameters
        + adata : AnnData An AnnData object
        + fname_prefix : str
        Prefix of the exported files. If not empty and not ending with '/' or
        '_', a '_' will be appended. Full names will be
        <fname_prefix>matrix.mtx, <fname_prefix>genes.tsv,
        <fname_prefix>barcodes.tsv
        + var : list
        A list of column names to be exported to gene table
        + obs : list
        A list of column names to be exported to barcode/cell table
    """
    if fname_prefix and not (
        fname_prefix.endswith("/") or fname_prefix.endswith("_")
    ):
        fname_prefix = fname_prefix + "_"
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
            mat = sp.coo_matrix(adata.layers[use_layer])
        else:
            mat = sp.coo_matrix(adata.X)

    obs = list(set(obs) & set(adata.obs.columns))
    var = list(set(var) & set(var_source.columns))

    n_obs, n_var = mat.shape
    n_entry = len(mat.data)
    header = (
        "%%MatrixMarket matrix coordinate real general\n%\n{} {} {}\n".format(
            n_var, n_obs, n_entry
        )
    )
    df = pd.DataFrame(
        {"col": mat.col + 1, "row": mat.row + 1, "data": mat.data}
    )
    mtx_fname = fname_prefix + "matrix.mtx"
    gene_fname = fname_prefix + "genes.tsv"
    barcode_fname = fname_prefix + "barcodes.tsv"
    with open(mtx_fname, "a") as fh:
        fh.write(header)
        df.to_csv(fh, sep=" ", header=False, index=False)

    obs_df = adata.obs[obs].reset_index(level=0)
    obs_df.to_csv(barcode_fname, sep="\t", header=False, index=False)
    var_df = var_source[var].reset_index(level=0)
    if not var:
        var_df["gene"] = var_df["index"]
    var_df.to_csv(gene_fname, sep="\t", header=False, index=False)


def calculate_summary_stats(adata, obs, matrix="normalised"):
    print(
        f"Calculating summary stats for {matrix} matrix, cell groups defined"
        f" by {obs}"
    )
    for ob in obs:
        layer = None
        use_raw = False

        if matrix == "raw.x":
            use_raw = True

        elif matrix in adata.layers:
            layer = matrix

        genedf = sc.get.obs_df(
            adata,
            keys=[ob, *list(adata.var_names)],
            layer=layer,
            use_raw=use_raw,
        )
        grouped = genedf.groupby(ob)
        mean, median = grouped.mean(), grouped.median()
        adata.varm[f"mean_{matrix}_{ob}"] = mean.transpose()
        adata.varm[f"median_{matrix}_{ob}"] = median.transpose()


if __name__ == "__main__":
    import doctest

    sys.exit(doctest.testmod(verbose=True)[0])
