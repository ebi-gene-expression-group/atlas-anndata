"""
Provides version, author and exports
"""

from .funcs import (
    validate_anndata_with_config,
    initialise_bundle,
    make_bundle_from_anndata,
)

__all__ = [
    "validate_anndata_with_config",
    "initialise_bundle",
    "make_bundle_from_anndata",
]
