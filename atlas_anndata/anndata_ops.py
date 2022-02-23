"""
Operations applied to annData objects, mostly extracting information in various
ways.
"""

import scanpy as sc
import re
import numpy as np
from collections import Counter
import pandas as pd
import scanpy_scripts as ss
from .util import (
    obs_markers,
)


def update_anndata(adata, config, matrix_for_markers=None):

    # Record the config in the object
    adata.uns["scxa_config"] = config

    # Reset the var names the the specified gene ID field
    adata.var.set_index(config["gene_meta"]["id_field"], inplace=True)

    # Calcluate markers where necessary
    marker_groupings = [
        x["slot"] for x in config["cell_meta"]["entries"] if x["markers"]
    ]
    if len(marker_groupings) > 0:
        calculate_markers(
            adata=adata, config=config, matrix=matrix_for_markers
        )


def overwrite_obs_with_magetab(adata, config, magetab_dir):
    """Overwrite the cell metadata for the anndata object with curated MAGE-TAB"""

    # Read the MAGE-TAB data, and probably remote [Comment] etc

    # For droplet, merge in sample-wise metdata, ensure order is as current
    # obs, and replace the obs dataframe.

    # For non-droplet, just synch the data frame and replace the obs

    return adata


def derive_metadata(adata, config, kind=None):

    """Extract cell metadata, select fields and separate out run-wise info."""

    sample_metadata = None
    cell_specific_metadata = None

    cell_metadata = adata.obs.copy()
    sample_field = sample_field = config["cell_meta"].get(
        "sample_field", "sample"
    )

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

        if sample_field in adata.obs.columns:

            print(
                "... deriving runs using supplied sample ID column"
                f" {sample_field}"
            )

            runs = adata.obs[sample_field]

            barcodes = [
                re.findall("[ATGC]{8,}", cell_id)[0]
                for sample_id, cell_id in zip(
                    adata.obs[sample_field], adata.obs_names
                )
            ]

        else:
            print("... deriving runs and barcodes by parsing cell IDs")
            runs, barcodes = parse_cell_ids(adata)

        # If we had to derive run IDs from cell IDs, save them in the metadata

        if sample_field is None:
            cell_metadata[sample_field] = runs

        # Add derived barcodes as a new column
        sample_colno = list(cell_metadata.columns).index(sample_field)
        cell_metadata.insert(
            (sample_colno + 1), column="barcode", value=barcodes
        )

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
            [sample_field]
            + [
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


def get_markers_table(adata, marker_grouping):

    de_table = ss.lib._diffexp.extract_de_table(
        adata.uns[f"markers_{marker_grouping}"]
    )
    de_table = de_table.loc[de_table.genes.astype(str) != "nan", :]

    return de_table


def calculate_markers(adata, config, matrix="X"):

    marker_groupings = [
        x["slot"] for x in config["cell_meta"]["entries"] if x["markers"]
    ]

    # If we're going to be doing any marker calculation, make sure the
    # indicated matrix is suitable

    if any([not obs_markers(adata, mg) for mg in marker_groupings]):

        layer = None

        matrix_description = [
            x for x in config["matrices"]["entries"] if x["slot"] == matrix
        ][0]

        if (
            matrix_description["scaled"]
            or (not matrix_description["normalised"])
            or (not matrix_description["log_transformed"])
        ):
            errmsg = (
                "Some markers need calculation, but the matrix indicated"
                f" ({matrix}) is not annotated in the input configuration as"
                " normalised, log transformed and unscaled as we would need"
                " for that. Please update annotations and/or perform matrix"
                " transformations as required."
            )
            raise Exception(errmsg)
        elif matrix != "X":
            layer = matrix

        for mg in marker_groupings:
            if not obs_markers(adata, mg):
                print(
                    f"Marker statistics not currently available for {mg},"
                    " recalculating with Scanpy..."
                )
                sc.tl.rank_genes_groups(
                    adata,
                    mg,
                    method="wilcoxon",
                    layer=layer,
                    key_added="markers_" + mg,
                )

    else:
        print()
