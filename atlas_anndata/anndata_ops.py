"""
Modifications to be made to annData objects themselves prior to Atlas
ingestion. This might be changing the scale of expression values,
re-calculating markers etc.
"""

import scanpy as sc
import re
import numpy as np
from collections import Counter
import pandas as pd


def derive_metadata(adata, config, kind=None):

    """Extract cell metadata, select fields and separate out run-wise info."""

    sample_metadata = None
    cell_specific_metadata = None

    cell_metadata = adata.obs.copy()
    sample_col = sample_col = config["cell_meta"].get("sample_col", None)

    # By default print all obs columns, but that's probably not we want in most
    # cases because of mixture of data types there (from curation, QC,
    # clustering etc)

    if kind is None:
        obs_columns = list(adata.obs.columns)
    else:
        obs_columns = [
            slot_def["slot"]
            for slot_def in config["cell_meta"]["entries"]
            if slot_def["kind"] == kind
        ]
        cell_metadata = cell_metadata[obs_columns]

    if config["droplet"]:

        # If a sample column has been supplied, then we can split the obs frame by
        # that, and determine the sample-wide metadata. We can also remove this
        # sample name from the cell identifier (often a sample/ barcode composite)
        # to make things tidier.

        print("Extracting metadata from annData object...")

        if sample_col is not None:

            print(
                "... deriving runs using supplied sample ID column"
                f" {sample_col}"
            )

            runs = adata.obs[sample_col]

            barcodes = [
                re.sub(
                    f"(.*[^-_ =/])[-_ =/]?{sample_id}[-_ =/]?(.*)",
                    "\\1",
                    cell_id,
                )
                for sample_id, cell_id in zip(
                    adata.obs[sample_col], adata.obs_names
                )
            ]

        else:
            print("... deriving runs and barcodes by parsing cell IDs")
            runs, barcodes = parse_cell_ids(adata)

        # If we had to derive run IDs from cell IDs, save them in the metadata

        if sample_col is None:
            adata.obs["sample"] = runs

        # Split cell metadata by run ID and create run-wise metadata with
        # any invariant value across all cells of a run

        print("... extracting metadata consistent within samples")
        unique_runs = list(set(runs))
        submetas = [cell_metadata[[y == x for y in runs]] for x in unique_runs]
        sample_metadata = pd.concat(
            [
                df[[x for x in df.columns if len(df[x].unique()) == 1]].head(1)
                for df in submetas
            ],
            join="inner",
        )
        sample_metadata["run"] = unique_runs
        sample_metadata.set_index("run", inplace=True)

        print("... assigning other metadata as cell_specific")
        cell_specific_metadata = cell_metadata[
            [
                x
                for x in cell_metadata.columns
                if x not in list(sample_metadata.columns)
            ]
        ]

    return cell_metadata, sample_metadata, cell_specific_metadata


def parse_cell_ids(adata, sample_name_col=None):
    """
    Cell names in droplet data are normally some composite of
    run/sample/library and barcode. But how that composite is done is the wild
    west. Maybe we can tame the madness.
    """

    id_patterns = {
        "atlas_standard": {
            "pattern": re.compile(r"(^\S+)-([ATGC]{8,})$"),
            "rearr_func": lambda string, regex: list(
                re.findall(regex, string)[0]
            ),
        },
        "just_barcode": {
            "pattern": re.compile(r"^([ATGC]+)$"),
            "rearr_func": lambda string, regex: re.findall(regex, string) * 2,
        },
        "barcode_at_end_with_sep": {
            "pattern": re.compile(r"(^\S+)[-_ =/]([ATGC]{8,})$"),
            "rearr_func": lambda string, regex: list(
                re.findall(regex, string)[0]
            ),
        },
        "barcode_at_start_with_sep": {
            "pattern": re.compile(r"^([ATGC]{8,})[-_ =/](\S+)$"),
            "rearr_func": lambda string, regex: list(
                re.findall(regex, string)[0]
            )[::-1],
        },
        "barcode_at_end_with_sep_and_suffix": {
            "pattern": re.compile(r"^(\S+)[-_ =/]([ATGC]{8,}\S+)$"),
            "rearr_func": lambda string, regex: list(
                re.findall(regex, string)[0]
            ),
        },
        "barcode_at_end_no_sep_with_suffix": {
            "pattern": re.compile(r"^(\S*[^ATGC])([ATGC]{8,}\S+)$"),
            "rearr_func": lambda string, regex: list(
                re.findall(regex, string)[0]
            ),
        },
    }

    def fix_id(cellid, pattern_type):
        fixer = id_patterns[pattern_type]["rearr_func"]
        pattern = id_patterns[pattern_type]["pattern"]
        return fixer(cellid, pattern)

    # Create a list with cell ID types according to some conventions

    cellid_types = np.array(["None"] * adata.n_obs, dtype=object)

    for pattern_name, pattern_settings in id_patterns.items():
        matching_names = [
            i
            for i, item in enumerate(adata.obs_names)
            if cellid_types[i] == "None"
            and re.search(pattern_settings["pattern"], item)
        ]
        cellid_types[matching_names] = pattern_name

    # Count how many IDs match each type. If we still have unknowns we'll need
    # a new pattern

    counts = dict(Counter(cellid_types))
    if "None" in counts:
        example_unknown = np.where(cellid_types == "None")[0][0]
        errmsg = (
            "Cannot identify all cell IDs as barcode/sample composities, e.g."
            f" {adata.obs_names[example_unknown]} for {example_unknown}th"
            " cell."
        )
        raise Exception(errmsg)

    else:

        id_parts = [
            fix_id(cellid, cellid_type)
            for cellid, cellid_type in zip(adata.obs_names, cellid_types)
        ]
        transposed = list(map(list, zip(*id_parts)))

        return transposed[0], transposed[1]
