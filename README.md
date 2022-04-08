# Single-cell Expression Atlas - compatible analysis bundles from annData files

In the Gene Expression team at the EBI we produce [Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/home) (SCXA), using a consistent pipeline to analyse data from the raw FASTQ files and produce the results visibile in the SCXA interface. Intermediate in this process is an 'analysis bundle' whereby the subset of results needed are gathered and formatted correctly for loading into our databases and indices. 

We sometimes encounter datasets without raw data, but where the user has an [annData](https://anndata.readthedocs.io/en/latest/) file containing expression matrices and derived results. The purpose of this repository is to provide a way of producing an analysis bundle directly from an annData object, allowing us (with some manual curation) to input these experiments. It will also be useful to simplify our own processes, since we can take the annData files now produced at the end of our analysis to produce a bundle in one step.

This is the implementation of an [internal strategy document](https://docs.google.com/document/d/1sdy9iOHKXUz8dEv66v1Cu77626_w2KuJ2myapWIN5So/edit#heading=h.mdsu1vbn6spl).

See also a walk-through for a specific [example dataset](EXAMPLE.md)

## Analysis bundles

### What make an analysis bundle?

An SCXA analysis bundle contains:

 - MTX-format expression matrices (raw, filtered, normalised)
 - TSV-format cell metadata used in analysis
 - TSV-format gene metadata used in analysis
 - TSV-format dimension reductions of different parameterisation (t-SNE, UMAP)
 - An annData-format file containing all the above
 - Reference files (GTF, cDNA) used in analysis
 - A software report detailing tools and versions used
 - A manifest file containing all the above

### MANIFEST

Here is an example manifest:

```
Description	File	Parameterisation
software_versions_file	software.tsv	
mtx_matrix_rows	filtered_normalised/genes.tsv.gz	filtered_normalised
mtx_matrix_cols	filtered_normalised/barcodes.tsv.gz	filtered_normalised
mtx_matrix_content	filtered_normalised/matrix.mtx.gz	filtered_normalised
tsv_matrix	filtered_normalised/filtered_normalised.tsv	filtered_normalised
mtx_matrix_rows	tpm/genes.tsv.gz	tpm
mtx_matrix_cols	tpm/barcodes.tsv.gz	tpm
mtx_matrix_content	tpm/matrix.mtx.gz	tpm
tsv_matrix	tpm/tpm.tsv	tpm
mtx_matrix_rows	raw/genes.tsv.gz	raw
mtx_matrix_cols	raw/barcodes.tsv.gz	raw
mtx_matrix_content	raw/matrix.mtx.gz	raw
tsv_matrix	raw/raw.tsv	raw
mtx_matrix_rows	raw_filtered/genes.tsv.gz	raw_filtered
mtx_matrix_cols	raw_filtered/barcodes.tsv.gz	raw_filtered
mtx_matrix_content	raw_filtered/matrix.mtx.gz	raw_filtered
tsv_matrix	raw_filtered/raw_filtered.tsv	raw_filtered
mtx_matrix_rows	tpm_filtered/genes.tsv.gz	tpm_filtered
mtx_matrix_cols	tpm_filtered/barcodes.tsv.gz	tpm_filtered
mtx_matrix_content	tpm_filtered/matrix.mtx.gz	tpm_filtered
tsv_matrix	tpm_filtered/tpm_filtered.tsv	tpm_filtered
cell_metadata	E-MTAB-6077.cell_metadata.tsv	
condensed_sdrf	E-MTAB-6077.condensed-sdrf.tsv	
project_file	E-MTAB-6077.project.h5ad	
reference_transcriptome	reference/Danio_rerio.GRCz11.cdna.all.104.fa.gz	
reference_annotation	reference/Danio_rerio.GRCz11.104.gtf.gz	
gene_metadata	reference/gene_annotation.txt	
protocol smart-seq
tsne_embeddings	tsne_perplexity_1.tsv	1
tsne_embeddings	tsne_perplexity_10.tsv	10
umap_embeddings	umap_n_neighbors_10.tsv	10
umap_embeddings	umap_n_neighbors_100.tsv	100
cluster_markers	markers_2.tsv	2
cluster_markers	markers_23.tsv	23
cluster_markers	markers_32.tsv	32
cluster_markers	markers_42.tsv	42
cluster_markers	markers_49.tsv	49
cluster_markers	markers_5.tsv	5
marker_stats	filtered_normalised_stats.csv	filtered_normalised
marker_stats	tpm_filtered_stats.csv	tpm_filtered
cluster_memberships	clusters_for_bundle.txt
```

## Process outline

The steps required to produce a bundle  are:

 1. Process the annData file to determine what information is available. Store a summary of this information in a YAML-format configuration file. 
 2. Generate bundle files (1st time) based on the starting config files in 1. 
 3. Examine the cell metadata files and use that information to refine configuration related to cell metadata. This includes flagging fields that should be included in the pre-MAGE-TAB files, and for droplet experiments, finding the field that separates cells from different libraries. Also check the gene metadata at this stage, to check that Ensembl gene identifiers are available.
 4. Generate bundle files (2nd time) based on configuration refined in 3.. This will include pre-MAGE-TAB files suitable as a basis for curation.
 5. Undertake curation. This will include standard curation processes used for SCXA-analysed experiments, and additional steps to derive further information we need, such as the nature of any included matrices and reference transcriptome used. At this point the YAML-format configuration should be complete, and should be included in the scxa-metadata repo.
 6. Generate the final bundle suitable for loading into SCXA.


## Protocol: producing an analysis bundle from arbitrary annData files

### 0. Install this package

This repository contains a python packages which should be used to facilitate production of a bundle. We'll get it on PyPi/ Conda soon, for now install like:

```
git clone git@github.com:ebi-gene-expression-group/atlas-anndata.git
cd atlas-anndata
pip install .
```

### 1. Produce a YAML format annData description file

To produce a valid bundle from an anndata file, we need to describe that file, outlining which of the cell/ gene metadata columns, matrices,dimension reductions etc should be included. This is done via a YAML-format config file (see [example](atlas_anndata/example_config.yaml)).

A starting configuration can be produced directly from the annData file using the `make_starting_config_from_anndata`. Usage for this command is:

```
Usage: make_starting_config_from_anndata [OPTIONS] ANNDATA_FILE ANNDATA_CONFIG

  Make a starting config describing an anndata file, for use in making
  analysis bundles for input into Single Cell Expression Atlas (SCXA)

  anndata_file   - A file of the annData hdf5 specification, with all
                   necessaryinformation for SCXA.
  anndata_config - File path to write YAML config.

Options:
  --atlas-style                  Assume the tight conventions from SCXA, e.g.
                                 on .obsm slot naming?

  --analysis-versions-file PATH  A four-column tab-delimited file with analys,
                                 analysis, version and citation

  --droplet                      Is this a droplet experiment?
  --gene-id-field TEXT           Field in .var where gene ID is stored.
  --gene-name-field TEXT         Field in .var where gene name (symbol) is
                                 stored.

  --sample-field TEXT            Field in .obs which separates cells from
                                 different libraries.

  --default-clustering TEXT      Of the unsupervised clusterings, which
                                 clustering should be set as the default? If
                                 not set, the middle (or first middle)
                                 clustering will be selected, or if --atlas-
                                 style is set, this will be the clustering
                                 corresponding to a resolution of 1.

  --help                         Show this message and exit.
```

(Note that the `atlas-style` flag is probably only useful for annData files produced by the Experession Atlas team, and relies on a number of assumptions about the content of the file in order to infer some additional information.)

For example to make a starting bundle for an annData file from a droplet experiment we might do:

```
make_starting_config_from_anndata project.h5ad test_config_from_anndata.yaml --droplet
``` 

This config is likely to wrong in a number of ways, but its just a starting point. 

### 2. Run bundle creation based on the unedited configuration file

Without editing the config we can do a first naive bundling run to output all the info we have. This is helpful in working out what information to complete in the config YAML file.

Detailed help for the relevant command is:

```
Usage: make_bundle_from_anndata [OPTIONS] ANNDATA_FILE ANNDATA_CONFIG
                                BUNDLE_DIR

  Build a bundle directory compatible with Single Cell Expression Atlas
  (SCXA) build proceseses

  anndata_file   - A file of the annData hdf5 specification, with all
                   necessaryinformation for SCXA.
  anndata_config - A config file generated with
                   `make_starting_config_from_anndata` and manually edited to
                   supply necessary information.
  bundle_dir     - A directory in which to create the bundle.

Options:
  --max-rank-for-stats INTEGER  For how many top marker genes should stats
                                (mean, median expression) be output?

  --matrix-for-markers TEXT     Where cell groups in the configuration file
                                have been flagged with markers, which matrix
                                should be used? Can be X, or an entry in
                                .layers(). The matrix must be appropriate for
                                Scanpy's tl.rank_genes_groups() method,
                                usually meaning filtered, normalised and log
                                transformed, but without additional scaling.
                                [required]

  --no-write-matrices           Use to disable writing of matrices to bundle,
                                for example in early metadata evaluation

  --write-premagetab            Should we write pre-magetab files for curation
                                by SCXA team? If not, we assume curation is
                                done, and we will look for curated metadata
                                via the the exp_name parameter

  --conda-prefix PATH           Specify a Conda directory to be used for
                                environments when running Snakemake workflows.

  --scxa-metadata-branch TEXT   When searching the SCXA metadata repository
                                for curation for this experiment, which branch
                                should we use?  [required]

  --sanitize-columns            When adding data from curation into the
                                anndata object, should we remove the Comment,
                                Characteristic etc?

  --exp-name TEXT               Specify an Expression Atlas identifier that
                                will be used for this experiment. If not set,
                                a placeholder value E-EXP-1234 will be used
                                and can be edited in the bundle later.
                                [required]

  --scxa-db-scale INTEGER       To what overall scale should cell counts be
                                multiplied for the SCXA DB? A multiplier will
                                be calculated from this value and the median
                                cell-wise sum in the given matrix.

  --help                        Show this message and exit.
```

For our first naive run we can supply basic arguments like:

```
make_bundle_from_anndata project.h5ad test_config_from_anndata.yaml test_bundle --write-premagetab
```

The `--write-premagetab` tells the package that this experiment is uncurated, so output some starting info for curation (though we won't use that yet).

### 3: refine configuration for curation

With the configuration file and unmoderated content available to us we can make some sensible decisions about some of those settings in the YAML file.

#### Identify correct gene ID and gene name fields

For SCXA we need the gene symbol and ID fields. The configuration YAML might have populated these if the default field names are present, but you may well get:

```
gene_meta:
  id_field: FILL ME with a string
  name_field: FILL ME with a string
```

You need to look at the `reference/gene_annotation.txt` file from the bundle directory and set these fields. `id_field` **must** be a field containing Ensembl gene IDs. If these are not available we cannot work with a dataset. `name_field` is a field containing gene symbols.

### Sample field

The sample field is encoded in the config generated above like:

```
  sample_field: FILL ME with a string
```

The value of this configuration field must be a field name from NONAME.cell_metadata.tsv in the bundle directory corresponding to a field that separates cells from different libraries.The `sample` field in this file is usually derived from the cell identifiers, and should only be used in the absence of more concrete information.

### Flag curated fields

Cell meta data from annData objects is a mixtrue of any input sample metadata provided by the author, plus annotations added over the course of analysis. The latter may not be appropriate for inclusion in the metadata in SCXA. Check the fields described in `cell_meta`, especially their kind ('curation', 'clustering', 'analysis'). Curated fields are those present before analysis, biological and technical info for cells and samples. Clustering is used to indicate the results of unsupervised cell clustering stored in .obs. Analysis is everthing else, comprising all other fields added to .obs during analysis.

Most importantly: 

 - Flip `curation` to `analysis` for any field entry which should not ultimately form part of the SCXA experiment MAGE-TAB format metadata.
 - Ensure any fields corresponding to unsupervised clusterings are flagged correctly.

### 4: second bundle run: re-write bundle starting metadata for curation

With the updated configuration in hand we can generate a bundle with pre-MAGE-TAB data for curation:

```
make_bundle_from_anndata project.h5ad test_config_from_anndata.yaml test_bundle --write-premagetab
```

This will produce the pre-MAGE-TAB files at `/mage-tab/NONAME.presdrf.txt` in the bundle.

### 5: Do curation

#### Regular curation efforts

The pre-MAGE-TAB can now be used to start curation by the curation team.

#### Gather other missing info needed before final bundle creation

Unlike our standard submission pathways, for pre-analysed data we need additional information before the data are ingested for SCXA, which must currently be provided via the configuration YAML. The completed config from the following steps should be added to the `scxa-metadata` alongside the MAGE-TAB files.

 - Under analyses please describe the analysis that was done. At a minimum you should describe the reference used (see [the example](atlas_anndata/example_config.yaml)) and the mapping tool used.
 - Under `matrices` check that you want all these matrices to be considered. You can remove any matrix that's not useful, and you should check the processing flags / matrices for each one. 
 - Under load_to_scxa_db please state the matrix that should be used by Atlas in expression-based displays. This should be filtered and normalised but not scaled or transformed. If no matrix in the object matches these criteria please remove this part of the config and Atlas will not show displays for this experiment based on expression values.
 - Check the dimension reductions described, again paying attention to 'kind'. 

For all sections, check [the example](atlas_anndata/example_config.yaml) for an idea of how things should look.For example under `matrices` there will be an entry pertaining the content of .X. You shold add a name (e.g 'scaled'), and c 

```
  - cell_filtered: true
    gene_filtered: true
    log_transformed: true
    measure: counts
    name: FILL ME with a string
    normalised: true
    parameters: {}
    scaled: true
    slot: X
```

You would fill the 'name' field here with something more descriptive for the matrix.

#### Validate config YAML

Having edited the config YAML, you should validate it against a schema we provide and the annData file itself. We can use this mechanism to ensure that inputs match the expectations of Single Cell Expression Atlas. 

```
validate_anndata_with_config test_config_from_anndata.yaml project.h5ad
```

There are no additional arguments for this step.

Bundling steps will also run this automatically before proceeding, but running it yourself will flag any issues early. If the validation flags any issues, resolve them.


### 6: third bundle run: write final bundle from completed YAML file

```
make_bundle_from_anndata --exp-name <exp name> --conda-prefix <conda prefix> project.h5ad test_config_from_anndata.yaml test_bundle
```

 - `exp name` is the E- accession generated for the experiment in curation
 - `conda prefix` is location conda environments will be stored for the workflow that does SDRF condensation etc.

This will pull the curated metadata from the `scxa-metadata` repo, condense the SDRF (adding ontology terms) and re-generate cell-wise annotations that will be used to enrich the content of the annData file. Other content such as matrices, dimensions etc will also be written, and this should form the final bundle that can be read by SCXA loading processes.

