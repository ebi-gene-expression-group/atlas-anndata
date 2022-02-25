import pkg_resources

MISSING = "FILL ME"
MISSING_STRING = f"{MISSING} with a string"
MISSING_BOOLEAN = f"{MISSING} with a boolean"

schema_file = pkg_resources.resource_filename(
    "atlas_anndata", "config_schema.yaml"
)
example_config_file = pkg_resources.resource_filename(
    "atlas_anndata", "example_config.yaml"
)

scxa_h5ad_test = pkg_resources.resource_filename(
    "atlas_anndata", "data/E-MTAB-6077.project.h5ad"
)


