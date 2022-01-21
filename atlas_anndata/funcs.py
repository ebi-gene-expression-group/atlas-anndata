import pkg_resources
import jsonschema
from jsonschema import validate
import yaml
import sys
import re
import scanpy as sc
import pandas as pd

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
            "leiden": "analysis",
            "log1p": "analysis",
            "^n_": "analysis",
            "total_": "analysis",
            "pct_": "analysis",
            "mito": "analysis",
            "louvain": "analysis",
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
            "name": ms if ms in adata.layers else None,
            "measure": "counts" if atlas_style else None,
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
        }

        markers_slot = obs_markers(adata, obs)
        if markers_slot:
            obs_entry["markers"] = True
            obs_entry["markers_slot"] = markers_slot
        else:
            obs_entry["markers"] = False

        config["cell_groups"]["entries"].append(obs_entry)

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
    gene_name_field=None,
):
    # Make sure the config matches the schema and anndata

    config, adata = validate_anndata_with_config(
        anndata_description_yaml, anndata_file
    )


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

    manifest = _read_file_manifest(bundle_dir)

    if write_cellmeta:
        print("Writing cell and run/sample metadata")
        manifest = _write_cell_metadata(
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
        manifest = _set_manifest_value(
            manifest, "gene_metadata", "reference/gene_annotation.txt"
        )

    if raw_matrix_slot is not None:
        print("Writing raw matrix")
        manifest = _write_matrix_from_adata(
            manifest,
            adata,
            slot=raw_matrix_slot,
            bundle_dir=bundle_dir,
            subdir="raw",
            gene_name_field=gene_name_field,
        )

    if filtered_matrix_slot is not None:
        print("Writing filtered matrix")
        manifest = _write_matrix_from_adata(
            manifest,
            adata,
            slot=filtered_matrix_slot,
            bundle_dir=bundle_dir,
            subdir="raw_filtered",
            gene_name_field=gene_name_field,
        )

    if normalised_matrix_slot is not None:
        print("Writing normalised matrix")
        manifest = _write_matrix_from_adata(
            manifest,
            adata,
            slot=normalised_matrix_slot,
            bundle_dir=bundle_dir,
            subdir="filtered_normalised",
            gene_name_field=gene_name_field,
        )

    if final_matrix_slot is not None:
        print("Writing final matrix")
        manifest = _write_matrix_from_adata(
            manifest,
            adata,
            slot=final_matrix_slot,
            bundle_dir=bundle_dir,
            subdir="final",
            gene_name_field=gene_name_field,
        )

    if write_obsms:
        print("Writing obsms to file")
        obsms = adata.obsm.keys() if obsms is None else obsms

        for obsm in obsms:
            print(f"Writing {obsm}")
            manifest = _write_obsm_from_adata(
                manifest,
                adata,
                obsm_name=obsm,
                bundle_dir=bundle_dir,
                atlas_style=atlas_style,
            )

    if write_clusters:
        print("Writing clusterings to file")

        # If no list of cluster variables is specified, select automatically
        # using the field pattern. If the fields are supplied, use the ones
        # that are actually present in .obs.

        if clusters is None:
            clusters = [
                x for x in adata.obs.columns if clusters_field_pattern in x
            ]
        else:
            clusters = list(set(clusters) & set(adata.obs.columns))

        if len(clusters) == 0:
            print("No .obs slots matching clusters specifications.")
            return False

        manifest = _write_clusters_from_adata(
            manifest,
            adata,
            clusters=clusters,
            bundle_dir=bundle_dir,
            atlas_style=atlas_style,
        )

    if write_markers:
        print("Writing markers to file")
        manifest = _write_markers_from_adata(
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
        manifest = _set_manifest_value(
            manifest, "software_versions_file", "software.tsv"
        )

    # Write anndata

    if write_anndata:

        print("Writing annData file to bundle")
        adata_filename = f"{config['exp-name']}.project.h5ad"
        adata.write(f"{bundle_dir}/{adata_filename}")
        manifest = _set_manifest_value(
            manifest, "project_file", adata_filename
        )

    # Write the final file manifest

    _write_file_manifest(bundle_dir, manifest)


# Write cell metadata, including for curation as mage-tab


def _write_cell_metadata(
    manifest, adata, bundle_dir, config, nonmeta_obs_patterns=None
):

    print("Writing cell metdata to be used in curation")
    cellmeta_filename = f"{bundle_dir}/{config['exp_name']}.cell_metadata.tsv"
    presdrf_filename = (
        f"{bundle_dir}/mage-tab/{config['exp_name']}.presdrf.txt"
    )
    precells_filename = (
        f"{bundle_dir}/mage-tab/{config['exp_name']}.precells.txt"
    )

    # Don't write calculated fields

    if nonmeta_obs_patterns is None:
        nonmeta_obs_patterns = [
            "louvain",
            "n_genes",
            "n_counts",
            "pct_",
            "total_counts",
            "predicted_doublet",
            "doublet_score",
        ]

    nonmeta_cols = [
        x
        for x in adata.obs.columns
        if any([y in x for y in nonmeta_obs_patterns])
    ]
    cellmeta_cols = [x for x in adata.obs.columns if x not in nonmeta_cols]

    cell_metadata = adata.obs[cellmeta_cols].copy()
    cell_metadata.to_csv(
        cellmeta_filename, sep="\t", header=True, index=True, index_label="id"
    )
    manifest = _set_manifest_value(
        manifest, "cell_metadata", cellmeta_filename
    )

    if config["droplet"]:

        # Split cell IDs to runs and barcodes

        runs, barcodes = zip(*(s.split("-") for s in adata.obs_names))
        cell_metadata["cell barcode"] = barcodes

        # Split cell metadata by run ID and create run-wise metadata with
        # any invariant value across all cells of a run

        unique_runs = list(set(runs))
        submetas = [cell_metadata[[y == x for y in runs]] for x in unique_runs]
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
            presdrf_filename,
            sep="\t",
            header=True,
            index=True,
            index_label="id",
        )
        cell_specific_metadata = cell_metadata[
            [
                x
                for x in cell_metadata.columns
                if x not in nonmeta_cols + list(run_meta.columns)
            ]
        ]

        if len(cell_specific_metadata.columns) > 0:
            cell_specific_metadata.to_csv(
                precells_filename,
                sep="\t",
                header=True,
                index=True,
                index_label="Cell ID",
            )
        else:
            print(
                "Supplied anndata contained no cell-specific metadata for this"
                " droplet experiment"
            )

    else:
        cell_metadata.to_csv(
            presdrf_filename,
            sep="\t",
            header=True,
            index=True,
            index_label="id",
        )

    return manifest


# Write dimension reductions from .obsm slots


def _write_obsm_from_adata(
    manifest, adata, obsm_name, bundle_dir, atlas_style=True
):

    obsm_paramstring = obsm_name.replace("X_", "")
    filename = f"{obsm_paramstring}.tsv"

    description = "embeddings"
    for embedding_type in ["umap", "pca", "tsne"]:
        if embedding_type in obsm_paramstring:
            description = f"{embedding_type}_{description}"
            obsm_paramstring = re.sub(
                r"_?%s_?" % embedding_type, "", obsm_paramstring
            )
            break

    # If following Atlas conventions, we can assume .obsm slots like
    # X_umap_neighbors_n_neighbors_20, which after above processing should now
    # be like neighbors_n_neighbors_20

    if atlas_style:
        parameterisation = obsm_paramstring.split("_")[-1]
    else:
        parameterisation = obsm_paramstring

    ss.obj_utils.write_embedding(
        adata, key=obsm_name, embed_fn=f"{bundle_dir}/{filename}"
    )

    manifest = _set_manifest_value(
        manifest, description, filename, parameterisation
    )

    return manifest


def _select_clusterings(adata, clusters, atlas_style=True):

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

    clustering_to_nclust = dict(
        zip(clusterings, [len(adata.obs[c].unique()) for c in clusterings])
    )

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


def _write_clusters_from_adata(
    manifest,
    adata,
    clusters,
    bundle_dir,
    atlas_style=True,
    default_clustering=None,
):

    if default_clustering is None:
        if atlas_style:
            default_clustering = "louvain_resolution_1.0"
        else:
            default_clustering = clusters[math.floor((len(clusters) - 1) / 2)]

    clustering_to_k = _select_clusterings(
        adata, clusters=clusters, atlas_style=atlas_style
    )

    # Write the complete clusters file in order

    with open(f"{bundle_dir}/clusters_for_bundle.txt", mode="w") as fp:
        fp.write("sel.K\tK\t" + "\t".join(list(adata.obs_names)) + "\n")

        for clustering, k in clustering_to_k.items():
            selected = "TRUE" if clustering == default_clustering else "FALSE"
            cluster_memberships = list(adata.obs[clustering])
            if min([int(x) for x in cluster_memberships]) == 0:
                cluster_memberships = [
                    str(int(x) + 1) for x in cluster_memberships
                ]

            fp.write("\t".join([selected, str(k)]))
            fp.write("\t" + "\t".join(list(cluster_memberships)))
            fp.write("\n")

    manifest = _set_manifest_value(
        manifest,
        "cluster_memberships",
        f"{bundle_dir}/clusters_for_bundle.txt",
    )

    return manifest


# For Atlas we store marker sets by integer number of clusters. Where cluster
# nubmers clash between marker sets (e.g. for different resolution values) we
# use the set from a resolution closest to 1
# TODO: fix this for cell type markers


def _write_markers_from_adata(
    manifest,
    adata,
    clusters,
    marker_clusterings=None,
    metadata_marker_fields=None,
    bundle_dir=None,
    atlas_style=True,
    write_marker_stats=True,
    max_rank_for_stats=4,
    marker_stats_layers=None,
):

    clustering_to_k = _select_clusterings(
        adata, clusters=clusters, atlas_style=atlas_style
    )

    # If no explicit marker sets suppplied, select those corresponding to the
    # selected clusterings (where available)

    marker_groupings = []

    if marker_clusterings is None:
        marker_groupings = [
            x
            for x in clustering_to_k.keys()
            if f"markers_{x}" in adata.uns.keys()
        ]

    if metadata_marker_fields is not None:
        marker_groupings = marker_groupings + metadata_marker_fields

    missing_marker_sets = [
        x for x in marker_groupings if f"markers_{x}" not in adata.uns.keys()
    ]
    if len(missing_marker_sets) > 0:
        raise Exception(
            "Some supplied marker clusterings do not have marker results in"
            " .uns: %s"
            % ",".join(missing_marker_sets)
        )

    de_tbls = dict(
        zip(
            marker_groupings,
            [_get_markers_table(adata, mg) for mg in marker_groupings],
        )
    )

    for mg in marker_groupings:
        marker = f"markers_{mg}"

        de_tbl = de_tbls[mg]

        # Reset cluster numbering to be from 1 if required

        if de_tbl["cluster"].min() == "0":
            de_tbl["cluster"] = [int(x) + 1 for x in de_tbl["cluster"]]

        if mg in clustering_to_k:
            k = clustering_to_k[mg]
            filename = f"{bundle_dir}/markers_{k}.tsv"
            de_tbl.to_csv(filename, sep="\t", header=True, index=False)
            manifest = _set_manifest_value(
                manifest, "cluster_markers", filename, k
            )
        else:
            filename = f"{bundle_dir}/{marker}.tsv"
            de_tbl.to_csv(filename, sep="\t", header=True, index=False)
            manifest = _set_manifest_value(
                manifest, "meta_markers", filename, mg
            )

    # Now make the summary stats

    if write_marker_stats:

        marker_stats_layers = (
            ["normalised"]
            if marker_stats_layers is None
            else marker_stats_layers
        )

        for sl in marker_stats_layers:

            if sl not in adata.layers.keys():
                valid_layers = ",".join(adata.layers.keys())
                raise Exception(
                    f"{sl} is not a valid layer, valid layers are"
                    f" {valid_layers}"
                )

            calculate_summary_stats(adata, marker_groupings, layer=sl)

            # Convert grouping name to k for storage
            marker_summary = pd.concat(
                [
                    _make_markers_summary(
                        adata,
                        sl,
                        mg,
                        de_tbl,
                        max_rank=max_rank_for_stats,
                        k=clustering_to_k.get(mg),
                    )
                    for mg, de_tbl in de_tbls.items()
                ]
            )
            statsfile = f"{bundle_dir}/{sl}_stats.csv"

            marker_summary.to_csv(statsfile, index=False)
            manifest = _set_manifest_value(
                manifest, "marker_stats", statsfile, sl
            )

    return manifest


def _get_markers_table(adata, marker_grouping):

    de_tbl = ss.lib._diffexp.extract_de_table(
        adata.uns[f"markers_{marker_grouping}"]
    )
    de_tbl = de_tbl.loc[de_tbl.genes.astype(str) != "nan", :]

    return de_tbl


def _make_markers_summary(
    adata, layer, marker_grouping, de_tbl, max_rank=4, k=None
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
        de_tbl.merge(summary_stats, left_on="genes", right_index=True)
        .drop(["ref", "scores", "logfoldchanges", "pvals"], axis=1)
        .rename(columns=new_colnames)
    )

    if max_rank:
        markers_summary = markers_summary[markers_summary["rank"] <= max_rank]

    # For unsupervised clusterings, record the grouping as k and increment the
    # group numbers so they start from 1

    if k:
        markers_summary["grouping_where_marker"] = k
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


def _read_file_manifest(bundle_dir):
    manifest_file = f"{bundle_dir}/MANIFEST"
    manifest = OrderedDict()

    if os.path.isfile(manifest_file):
        with open(manifest_file) as fp:
            _ = fp.readline()
            for line in fp:
                line_parts = line.rstrip().split("\t")
                if len(line_parts) < 3:
                    line_parts.append("")
                manifest = _set_manifest_value(
                    manifest, line_parts[0], line_parts[1], line_parts[2]
                )

    return manifest


def _write_file_manifest(bundle_dir, manifest):
    manifest_file = f"{bundle_dir}/MANIFEST"

    with open(manifest_file, "w") as fh:
        fh.write("Description\tFile\tParameterisation\n")

        for description, v in manifest.items():
            for parameterisation, filename in v.items():
                fh.write(f"{description}\t{filename}\t{parameterisation}\n")


def _set_manifest_value(manifest, description, filename, parameterisation=""):
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


def _write_matrix_from_adata(
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

    manifest = _set_manifest_value(
        manifest, "mtx_matrix_content", f"{subdir}/matrix.mtx.gz", subdir
    )
    manifest = _set_manifest_value(
        manifest, "mtx_matrix_cols", f"{subdir}/barcodes.tsv.gz", subdir
    )
    manifest = _set_manifest_value(
        manifest, "mtx_matrix_rows", f"{subdir}/genes.mtx.gz", subdir
    )

    return manifest


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


def calculate_summary_stats(adata, obs, layer="normalised"):
    print("Calculating summary stats")
    for ob in obs:
        genedf = sc.get.obs_df(
            adata, keys=[ob, *list(adata.var_names)], layer=layer
        )
        grouped = genedf.groupby(ob)
        mean, median = grouped.mean(), grouped.median()
        adata.varm[f"mean_{layer}_{ob}"] = mean.transpose()
        adata.varm[f"median_{layer}_{ob}"] = median.transpose()


if __name__ == "__main__":
    import doctest

    sys.exit(doctest.testmod(verbose=True)[0])
