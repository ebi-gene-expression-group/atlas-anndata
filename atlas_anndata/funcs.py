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
import csv
from .anndata_ops import derive_metadata

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

    # Validate against the schema
    schema = load_doc(schema_file)

    try:
        validate(instance=config, schema=schema)
    except jsonschema.exceptions.ValidationError as err:
        print(err, file=sys.stderr)
        return False

    # Also check that no blank values have been left in

    from atlas_anndata import MISSING

    if MISSING in str(config):
        errmsg = (
            f"Please complete all {MISSING} fields in config before trying to"
            " make a bundle."
        )
        raise Exception(errmsg)

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

    check_result = False

    errmsg = f"{slot_type} entry {slot_name} not present in anndata file"

    if slot_type == "matrices":
        if slot_name == "X":
            check_result = hasattr(adata, "X")
        elif slot_name == "raw.X":
            check_result = hasattr(adata.raw, "X")
        else:
            check_result = slot_name in adata.layers

    elif slot_type == "dimension_reductions":
        check_result = slot_name in adata.obsm

    elif slot_type == "cell_meta":
        check_result = slot_name in adata.obs.columns or slot_name == "index"

    elif slot_type == "gene_meta":
        check_result = slot_name in adata.var.columns or slot_name == "index"

    else:
        errmsg = f"{slot_type} slot type not recognised"
        check_result = False

    if not check_result:
        raise Exception(errmsg)

    return check_result


def validate_anndata_with_config(anndata_config, anndata_file):

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
    Checking for cell_meta organism_part
    Checking for cell_meta louvain_resolution_0.7
    Checking for cell_meta louvain_resolution_1.0
    Checking for dimension_reductions X_umap_neighbors_n_neighbors_3
    Checking for dimension_reductions X_umap_neighbors_n_neighbors_10
    Checking for dimension_reductions X_umap_neighbors_n_neighbors_10
    Checking for gene_meta gene_name
    annData file successfully validated against config ...
    """

    config = load_doc(anndata_config)

    # First validate the anndata descripton file against the YAML schema

    print(f"Validating {anndata_config} against {schema_file}")
    config_status = validate_config(config)

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
                check_slot(adata, slot_type, slot_def["slot"])

    # Check that scxa load matrix (if specified) is present

    if "load_to_scxa_db" in config["matrices"]:
        check_slot(adata, "matrices", config["matrices"]["load_to_scxa_db"])

    check_slot(adata, "gene_meta", config["gene_meta"]["id_field"])
    check_slot(adata, "gene_meta", config["gene_meta"]["name_field"])

    if "sample_field" in config["cell_meta"]:
        check_slot(adata, "cell_meta", config["cell_meta"]["sample_field"])

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
        kind = f"{kind}: 'pca', 'tsne', 'umap', 'scanvi'"
        search_map = {
            ".*pca": "pca",
            ".*umap": "umap",
            ".*tsne": "tsne",
            ".*scanvi.*": "scanvi",
        }

    elif slot_type == "cell_meta":
        kind = "curation"
        search_map = {
            ".*leiden.*": "clustering",
            "cluster": "clustering",
            "log1p": "analysis",
            "^n_": "analysis",
            "total_": "analysis",
            "pct_": "analysis",
            "mito": "analysis",
            ".*louvain.*": "clustering",
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
    ... 'cell_meta',
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

        elif slot_type == "cell_meta":
            m = re.search(r".*(louvain|leiden)_(.*)_(.*)", slot_name)
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


def read_analysis_versions_file(analysis_versions_file, atlas_style=False):

    analysis_versions = pd.read_csv(analysis_versions_file, sep="\t")
    analysis_versions.columns = [x.lower() for x in analysis_versions.columns]

    # We need a 'kind' to differentiate software and e.g. reference files.
    # Atlas hasn't had it up to now (we kind of shoehorned references in
    # to the software field), but we know what the content of our versions
    # file is, so we can infer it.

    analysis_versions.rename(columns={"software": "asset"}, inplace=True)

    required_versions_columns = [
        "analysis",
        "asset",
        "version",
        "citation",
    ]
    if not atlas_style:
        required_versions_columns.append("kind")

    if not set(required_versions_columns).issubset(analysis_versions.columns):
        errmsg = (
            "Analysis versions file must have at least columns"
            f" {required_versions_columns}, recieved"
            f" {analysis_versions.columns}"
        )
        raise Exception(errmsg)
    else:
        if "kind" not in analysis_versions.columns:
            required_versions_columns.append("kind")
            analysis_versions["kind"] = [
                "file" if x.lower() == "reference" else "software"
                for x in analysis_versions["analysis"]
            ]

    # Remove any other columns that might be present

    analysis_versions = analysis_versions[required_versions_columns]

    return analysis_versions


def make_starting_config_from_anndata(
    anndata_file,
    anndata_config,
    atlas_style=False,
    exp_name=None,
    droplet=False,
    gene_id_field="gene_id",
    gene_name_field="gene_name",
    sample_field="sample",
    default_clustering=None,
    analysis_versions_file=None,
):

    """
    Make a yaml-format configuration file as a starting point for manual
    editing, from the content of a provided annData file.

    >>> make_starting_config_from_anndata(scxa_h5ad_test, '/tmp/foo.yaml')
    """

    adata = sc.read(anndata_file)
    from atlas_anndata import MISSING_STRING

    config = {
        "droplet": droplet,
        "matrices": {"load_to_scxa_db": MISSING_STRING, "entries": []},
        "gene_meta": {
            "name_field": gene_name_field
            if gene_name_field in adata.var.columns
            else MISSING_STRING,
            "id_field": gene_id_field
            if gene_id_field in adata.var.columns
            else MISSING_STRING,
        },
        "cell_meta": {"entries": []},
        "dimension_reductions": {"entries": []},
        "analysis_versions": [],
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

    # Check that we actually have some obs

    if len(adata.obs.columns) == 0:
        errmsg = (
            "Object in {anndata_file} has no obs (cell metadata at all) and as"
            " such is not a candidate for inclusion in SCXA."
        )
        raise Exception(errmsg)

    # Describe cell-wise metadata columns

    for obs in adata.obs.columns:

        obs_entry = {
            "slot": obs,
            "kind": slot_kind_from_name("cell_meta", obs),
            "parameters": extract_parameterisation(
                "cell_meta", obs, atlas_style
            ),
            "default": False,
        }

        if default_clustering is not None and obs == default_clustering:
            obs_entry["default"] = default_clustering

        markers_slot = obs_markers(adata, obs)
        if markers_slot:
            obs_entry["markers"] = True
            obs_entry["markers_slot"] = markers_slot
        else:
            obs_entry["markers"] = False

        config["cell_meta"]["entries"].append(obs_entry)

    # For Atlas objects we'll use the predictable cell ID structure to derive
    # samples and barcodes. For other sources we probably need an explicit
    # sample column

    if droplet and not atlas_style:
        config["cell_meta"]["sample_field"] = (
            sample_field
            if sample_field in adata.obs.columns
            else MISSING_STRING
        )

    # Find the groups in obs that correspond to clusterings
    config_obs = [x["slot"] for x in config["cell_meta"]["entries"]]
    cluster_obs = [
        x["slot"]
        for x in config["cell_meta"]["entries"]
        if x["kind"] == "clustering"
    ]

    cluster_obs = list(
        select_clusterings(
            adata, clusters=cluster_obs, atlas_style=atlas_style
        ).keys()
    )

    if default_clustering is None:

        # Try to set a default clustering
        if atlas_style and "louvain_resolution_1.0" in cluster_obs:
            default_clustering = "louvain_resolution_1.0"
        else:
            default_clustering = cluster_obs[
                math.floor((len(cluster_obs) - 1) / 2)
            ]

        config["cell_meta"]["entries"][config_obs.index(default_clustering)][
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

    # Add some placeholders to encourage users to fill in software

    if analysis_versions_file is not None:
        analysis_versions = read_analysis_versions_file(
            analysis_versions_file, atlas_style=atlas_style
        )
        config["analysis_versions"] = list(
            analysis_versions.T.to_dict().values()
        )
    else:
        for analysis in [
            "reference",
            "filtering and trimming",
            "mapping",
            "clustering",
        ]:
            config["analysis_versions"].append(
                {
                    "analysis": analysis,
                    "kind": "file" if analysis == "reference" else "software",
                    "asset": MISSING_STRING,
                    "version": MISSING_STRING,
                    "citation": MISSING_STRING,
                }
            )

    with open(anndata_config, "w") as file:
        yaml.dump(config, file)


def make_bundle_from_anndata(
    anndata_file,
    anndata_config,
    bundle_dir,
    max_rank_for_stats=5,
    exp_name="NONAME",
    write_premagetab=False,
    write_matrices=True,
    **kwargs,
):
    # Make sure the config matches the schema and anndata

    config, adata = validate_anndata_with_config(anndata_config, anndata_file)

    # Clear and create the output location

    if Path(bundle_dir).is_dir():
        shutil.rmtree(bundle_dir)

    pathlib.Path(f"{bundle_dir}").mkdir(parents=True)

    # Initialise the manifest

    manifest = read_file_manifest(bundle_dir)

    # Write cell metadata (curated cell info)

    print("Writing obs (curated metadata)")
    write_cell_metadata(
        manifest=manifest,
        adata=adata,
        bundle_dir=bundle_dir,
        config=config,
        kind="curation",
        exp_name=exp_name,
        write_premagetab=write_premagetab,
    )

    # Write matrices

    if write_matrices:
        print("Writing matrices")
        write_matrices_from_adata(
            manifest=manifest,
            bundle_dir=bundle_dir,
            adata=adata,
            config=config,
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

    # Write any associated markers

    print("Writing markers and statistics")
    write_markers_from_adata(
        manifest=manifest,
        bundle_dir=bundle_dir,
        adata=adata,
        config=config,
        write_marker_stats=True,
        max_rank_for_stats=max_rank_for_stats,
    )

    # Write any dim. reds from obsm

    print("Writing dimension reductions")
    write_obsms_from_adata(
        manifest=manifest,
        bundle_dir=bundle_dir,
        adata=adata,
        config=config,
    )

    # Write software table

    print("Writing analysis versions table")
    if len(config["analysis_versions"]) > 0:
        pd.DataFrame(config["analysis_versions"]).to_csv(
            f"{bundle_dir}/software.tsv", sep="\t", index=False
        )
        set_manifest_value(
            manifest, "analysis_versions_file", f"{bundle_dir}/software.tsv"
        )

    print("Writing annData file")

    # Record the config in the object
    adata.uns["scxa_config"] = config

    adata_filename = f"{exp_name}.project.h5ad"
    adata.write(f"{bundle_dir}/{adata_filename}")
    set_manifest_value(manifest, "project_file", adata_filename)

    # Write the final file manifest

    write_file_manifest(bundle_dir, manifest)


def write_matrices_from_adata(
    manifest,
    bundle_dir,
    adata,
    config,
):
    for slot_def in config["matrices"]["entries"]:
        write_matrix_from_adata(
            manifest=manifest,
            adata=adata,
            slot=slot_def["slot"],
            bundle_dir=bundle_dir,
            subdir=slot_def["name"],
            gene_name_field=config["gene_meta"]["name_field"],
        )


def write_clusters_from_adata(manifest, bundle_dir, adata, config):

    # Find the groups in obs that correspond to clusterings
    cluster_obs = [
        x["slot"]
        for x in config["cell_meta"]["entries"]
        if x["kind"] == "clustering"
    ]
    default_cluster_obs = [
        x["default"]
        for x in config["cell_meta"]["entries"]
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

    print("Writing cell metdata to be used in curation")
    cellmeta_filename = f"{exp_name}.cell_metadata.tsv"
    presdrf_filename = f"mage-tab/{exp_name}.presdrf.txt"
    precells_filename = f"mage-tab/{exp_name}.precells.txt"
    pathlib.Path(f"{bundle_dir}/mage-tab").mkdir(parents=True, exist_ok=True)

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
    max_rank_for_stats=5,
):
    marker_groupings_kinds = [
        (x["slot"], x["kind"])
        for x in config["cell_meta"]["entries"]
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
                for cell_grouping, de_table in de_tables.items()
            ]
        )
        statsfile = f"{bundle_dir}/{matrix_for_stats}_stats.csv"

        marker_summary.to_csv(
            statsfile, index=False, quoting=csv.QUOTE_NONNUMERIC
        )
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
    adata, layer, marker_grouping, de_table, max_rank=5, cell_group_kind=None
):

    print(f"... calculating stats for cell grouping {marker_grouping}")

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
        print(
            f"...... limiting stats report to top {max_rank} differential"
            " genes"
        )
        markers_summary = markers_summary[
            markers_summary["rank"] <= (max_rank - 1)
        ]

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
