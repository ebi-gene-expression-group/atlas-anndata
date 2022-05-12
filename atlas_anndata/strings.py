import pkg_resources

MISSING = "FILL ME"
MISSING_STRING = f"{MISSING} with a string"
MISSING_BOOLEAN = f"{MISSING} with a boolean"

load_to_dbxa_matrix_name = "filtered_normalised"

schema_file = pkg_resources.resource_filename(
    "atlas_anndata", "config_schema.yaml"
)
example_bundle_dir = pkg_resources.resource_filename(
    "atlas_anndata", "data/bundles"
)
example_config_file = pkg_resources.resource_filename(
    "atlas_anndata", "data/bundles/E-MTAB-6077/anndata-config.yaml"
)
example_software_file = pkg_resources.resource_filename(
    "atlas_anndata", "example_software.tsv"
)
example_manifest = pkg_resources.resource_filename(
    "atlas_anndata", "data/bundles/E-MTAB-6077/MANIFEST"
)
scxa_h5ad_test = pkg_resources.resource_filename(
    "atlas_anndata", "data/bundles/E-MTAB-6077/E-MTAB-6077.project.h5ad"
)
workflow_dir = pkg_resources.resource_filename("atlas_anndata", "snakemake")
