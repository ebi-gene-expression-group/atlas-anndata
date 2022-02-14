"""
Provides version, author and exports
"""

from .funcs import (
    validate_anndata_with_config,
    make_starting_config_from_anndata,
    make_bundle_from_anndata,
)

__all__ = [
    "validate_anndata_with_config",
    "make_starting_config_from_anndata",
    "make_bundle_from_anndata",
]

MISSING = "FILL ME"
MISSING_STRING = f"{MISSING} with a string"
MISSING_BOOLEAN = f"{MISSING} with a boolean"
