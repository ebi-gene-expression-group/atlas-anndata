import pkg_resources
import jsonschema
from jsonschema import validate
import yaml
import sys
import re
import scanpy as sc

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
            for slot_def in config[slot_type]:
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
        "matrices": [],
        "cell_groups": [],
        "dimension_reductions": list(),
    }

    # Describe the matrices

    matrix_slots = []

    if hasattr(adata, "raw") and hasattr(adata.raw, "X"):
        matrix_slots.append("raw.X")

    matrix_slots = matrix_slots + list(adata.layers.keys())
    matrix_slots.append("X")

    for ms in matrix_slots:
        matrix_entry = {
            "slot": ms,
            "measure": "counts",
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

        config["matrices"].append(matrix_entry)

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

        config["cell_groups"].append(obs_entry)

    # Describe dimension reductions stored in .obsm

    for obsm in adata.obsm.keys():
        obsm_entry = {
            "slot": obsm,
            "kind": slot_kind_from_name("dimension_reductions", obsm),
            "parameters": extract_parameterisation(
                "dimension_reductions", obsm, atlas_style
            ),
        }

        config["dimension_reductions"].append(obsm_entry)

    with open(config_file, "w") as file:
        yaml.dump(config, file)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
