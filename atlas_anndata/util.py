from .strings import MISSING_STRING
import scanpy as sc
import re
import yaml
import scanpy_scripts as ss
import pandas as pd
from .strings import (
    schema_file,
    example_config_file,
    example_software_file,
    scxa_h5ad_test,
)


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


def clusterings_to_ks(adata, obs_names):

    """Make a dictionary mapping specified cell groupings to numbers of clusters.

    >>> adata = sc.read(scxa_h5ad_test)
    >>> from atlas_anndata.anndata_config import load_doc
    >>> egconfig = load_doc(example_config_file)
    >>> cluster_obs = [
    ... x["slot"]
    ... for x in egconfig["cell_meta"]["entries"]
    ... if x["kind"] == "clustering"
    ... ]
    >>> clusterings_to_ks(adata, cluster_obs)
    {'louvain_resolution_0.7': 2, 'louvain_resolution_1.0': 5}
    """

    return dict(
        zip(obs_names, [len(adata.obs[c].unique()) for c in obs_names])
    )


def select_clusterings(adata, clusterings, atlas_style=True):

    """
    Pick clusterings to use. This is a wrapper around clusterings_to_ks(), for
    example to remove multiple clusterings with the same k value

    >>> adata = sc.read(scxa_h5ad_test)
    >>> from atlas_anndata.anndata_config import load_doc
    >>> egconfig = load_doc(example_config_file)
    >>> cluster_obs = [
    ... x["slot"]
    ... for x in egconfig["cell_meta"]["entries"]
    ... if x["kind"] == "clustering"
    ... ]
    >>> select_clusterings(adata, cluster_obs, atlas_style = True)
    {'louvain_resolution_0.7': 2, 'louvain_resolution_1.0': 5}
    """

    if atlas_style:
        clusterings = list(
            dict(
                sorted(
                    dict(
                        zip(
                            clusterings,
                            [
                                abs(1 - float(c.split("_")[-1]))
                                for c in clusterings
                            ],
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


def check_slot(adata, slot_type, slot_name):

    """Check for a slot in an anndata object

    >>> adata = sc.read(scxa_h5ad_test)
    >>> check = check_slot(adata, 'matrices', 'X')
    ..Checking for matrices X
    >>> check
    True
    """

    print(f"..Checking for {slot_type} {slot_name}")

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


def slot_kind_from_name(slot_type, slot_name):

    """
    Try to infer the kind/sub-type of a slot by comparing to known patterns


    >>> slot_kind_from_name('dimension_reductions', 'X_tsne_blah')
    'tsne'
    """

    kind = MISSING_STRING
    search_map = {}

    if slot_type == "dimension_reductions":
        kind = "other"
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


def read_analysis_versions_file(analysis_versions_file, atlas_style=False):

    """Read a tab-separated versions file into a Pandas dataframe

    >>> analysis = read_analysis_versions_file(example_software_file, atlas_style = True)
    >>> analysis.head()
                   analysis  ...      kind
    0            Build List  ...  software
    1  clustering_slotnames  ...  software
    2        Column arrange  ...  software
    3                   Cut  ...  software
    4          filter_cells  ...  software
    <BLANKLINE>
    [5 rows x 5 columns]
    """

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
