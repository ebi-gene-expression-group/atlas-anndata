import pkg_resources
import yaml
import jsonschema
from jsonschema import validate
from .strings import MISSING, MISSING_STRING
import sys

from .util import (
    slot_kind_from_name,
    extract_parameterisation,
    obs_markers,
    select_clusterings,
    read_analysis_versions_file,
)
import math

schema_file = pkg_resources.resource_filename(
    "atlas_anndata", "config_schema.yaml"
)
example_config_file = pkg_resources.resource_filename(
    "atlas_anndata", "example_config.yaml"
)


def describe_matrices(adata, atlas_style=False, droplet=False):

    conf = {"load_to_scxa_db": MISSING_STRING, "entries": []}

    # Describe the matrices

    matrix_slots = []

    if hasattr(adata, "raw") and hasattr(adata.raw, "X"):
        matrix_slots.append("raw.X")

    matrix_slots = matrix_slots + list(adata.layers.keys())
    matrix_slots.append("X")

    if atlas_style and "normalised" in matrix_slots:
        conf["load_to_scxa_db"] = "normalised"

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

        conf["entries"].append(matrix_entry)

    return conf


def describe_cellmeta(
    adata,
    atlas_style=False,
    droplet=False,
    default_clustering=None,
    sample_field="sample",
):

    conf = {"entries": []}

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

        conf["entries"].append(obs_entry)

    # For Atlas objects we'll use the predictable cell ID structure to derive
    # samples and barcodes. For other sources we probably need an explicit
    # sample column

    if droplet and not atlas_style:
        conf["sample_field"] = (
            sample_field
            if sample_field in adata.obs.columns
            else MISSING_STRING
        )

    # Find the groups in obs that correspond to clusterings
    config_obs = [x["slot"] for x in conf["entries"]]
    cluster_obs = [
        x["slot"] for x in conf["entries"] if x["kind"] == "clustering"
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

        conf["entries"][config_obs.index(default_clustering)]["default"] = True

    return conf


def describe_dimreds(adata, atlas_style=False, droplet=False):

    conf = {"entries": []}

    # Describe dimension reductions stored in .obsm

    for obsm in adata.obsm.keys():
        obsm_entry = {
            "slot": obsm,
            "kind": slot_kind_from_name("dimension_reductions", obsm),
            "parameters": extract_parameterisation(
                "dimension_reductions", obsm, atlas_style
            ),
        }

        conf["entries"].append(obsm_entry)

    return conf


def describe_genemeta(
    adata,
    atlas_style=False,
    droplet=False,
    gene_id_field="gene_id",
    gene_name_field="gene_name",
):
    return {
        "name_field": gene_name_field
        if gene_name_field in adata.var.columns
        else MISSING_STRING,
        "id_field": gene_id_field
        if gene_id_field in adata.var.columns
        else MISSING_STRING,
    }


def describe_analysis(
    adata, atlas_style=False, droplet=False, analysis_versions_file=None
):

    conf = []

    # Add some placeholders to encourage users to fill in software

    if analysis_versions_file is not None:
        analysis_versions = read_analysis_versions_file(
            analysis_versions_file, atlas_style=atlas_style
        )
        conf = list(analysis_versions.T.to_dict().values())
    else:
        for analysis in [
            "reference",
            "filtering and trimming",
            "mapping",
            "clustering",
        ]:
            conf.append(
                {
                    "analysis": analysis,
                    "kind": "file" if analysis == "reference" else "software",
                    "asset": MISSING_STRING,
                    "version": MISSING_STRING,
                    "citation": MISSING_STRING,
                }
            )

    return conf

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
    
    print(f"Validating config against {schema_file}")

    try:
        validate(instance=config, schema=schema)
    except jsonschema.exceptions.ValidationError as err:
        print(err, file=sys.stderr)
        return False

    # Also check that no blank values have been left in

    if MISSING in str(config):
        errmsg = (
            f"Please complete all {MISSING} fields in config before trying to"
            " make a bundle."
        )
        raise Exception(errmsg)

    return True
