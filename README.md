# Single-cell Expression Atlas - compatible analysis bundles from annData files

In the Gene Expression team at the EBI we produce [Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/home) (SCXA), using a consistent pipeline to analyse data from the raw FASTQ files and produce the results visibile in the SCXA interface. Intermediate in this process is an 'analysis bundle' whereby the subset of results needed are gathered and formatted correctly for loading into our databases and indices. 

We sometimes encounter datasets without raw data, but where the user has an [annData](https://anndata.readthedocs.io/en/latest/) file containing expression matrices and derived results. The purpose of this repository is to provide a way of producing an analysis bundle directly from an annData object, allowing us (with some manual curation) to input these experiments. It will also be useful to simplify our own processes, since we can take the annData files now produced at the end of our analysis to produce a bundle in one step.

This is the implementation of an [internal strategy document](https://docs.google.com/document/d/1sdy9iOHKXUz8dEv66v1Cu77626_w2KuJ2myapWIN5So/edit#heading=h.mdsu1vbn6spl).

## Analysis bundles

### What make analysis bundle?

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

## Protocol: producing an analysis bundle from arbitrary annData files

### 0. Install this package

This repository contains a python packages which should be used to facilitate production of a bundle. We'll get it on PyPi/ Conda soon, for now install like:

```
git clone git@github.com:ebi-gene-expression-group/atlas-anndata.git
cd atlas-anndata
pip install .
```

### 1. Produce a YAML format annData description file

To produce a valid bundle from an anndata file, we need to describe that file, outlining which of the cell/ gene metadata columns, matrices,dimension reductions etc should be included. This is done via a YAML-format config file (see [example]{example_config.yaml}).

A starting configuration can be produced directly from the annData file using the `make_starting_config_from_anndata`. Usage for this command is:

```
Usage: make_starting_config_from_anndata [OPTIONS] ANNDATA_FILE ANNDATA_CONFIG

  Build a bundle directory compatible with Single Cell Expression Atlas (SCXA)
  build proceseses

  anndata_file   - A file of the annData hdf5 specification, with all
                   necessaryinformation for SCXA.
  anndata_config - File path to write YAML config.

Options:
  --atlas-style                  Assume the tight conventions from SCXA, e.g.
                                 on .obsm slot naming?
  --software-versions-file PATH  A four-column tab-delimited file with analys,
                                 software, version and citation
  --droplet                      Is this a droplet experiment?
  --gene-name-field TEXT         Field in .var where gene name (symbol) is
                                 stored.
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

#### Complete the information in the YAML config file

The resulting YAML will be populated decriptions of the content in the various parts of a YAML file. This will largely be the result of guesswork, and you should check all content and update accordingly. Specifically:

 - Under analyses please desribe the analysis that was done. At a minimum you should describe the reference used (see [the example](atlas_anndata/example_config.yaml)) and the mapping tool used.
 - Under `matrices` check that you want all these matrices to be considered. You can remove any matrix that's not useful, and you should check the processing flags / matrices for each one. 
 - Under load_to_scxa_db please state the matrix that should be used by Atlas in expression-based displays. This should be filtered and normalised but not scaled or transformed. If no matrix in the object matches these criteria please remove this part of the config and Atlas will not show displays for this experiment based on expression values.
 - Under `gene_meta` state the field in .var which can be used as gene name/symbol. In Ensembl data this is conventionally 'gene_name'.
 - Check the fields described in `cell_meta`, especially their kind ('curation', 'clustering', 'analysis'). Curated fields are those present before analysis, biological and technical info for cells and samples. Clustering is used to indicate the results of unsupervised cell clustering stored in .obs. Analysis is everthing else, comprising all other fields added to .obs during analysis.
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

### 2. Validate config YAML

Having edited the config YAML, you should validate it against a schema we provide and the annData file itself. We can use this mechanism to ensure that inputs match the expectations of Single Cell Expression Atlas. 

```
validate_anndata_with_config test_config_from_anndata.yaml project.h5ad
```

There are no additional arguments for this step.

Subsequent steps will also run this automatically before proceeding, but running it yourself will flag any issues early. If the validation flags any issues, resolve them.

### 3. Run bundle creation based on the supplied configuration

Once the configuration file is complete, with all necesary info, the bundle generation itself can be done.

Detailed help for the relevant command is:

```
Usage: make_bundle_from_anndata [OPTIONS] ANNDATA_FILE ANNDATA_CONFIG
                                BUNDLE_DIR

  Build a bundle directory compatible with Single Cell Expression Atlas (SCXA)
  build proceseses

  anndata_file   - A file of the annData hdf5 specification, with all
                   necessaryinformation for SCXA.
  anndata_config - A config file generated with
                   `make_starting_config_from_anndata` and manually edited to
                   supply necessary information.
  bundle_dir     - A directory in which to create the bundle.

Options:
  --max-rank-for-stats INTEGER  For how many top marker genes should stats
                                (mean, median expression) be output?
  --write-premagetab            Should we write pre-magetab files for curation
                                by SCXA team?
  --magetab-dir PATH
  --exp-name TEXT               Specify an Expression Atlas identifier that
                                will be used for this experiment. 
  --help                        Show this message and exit.
```

### 3a. First run to derive starting metadata

The initial step for getting data represented in an annData object into SCXA is generation of full MAGE-TAB format metadata. We can generate starting content for the curators like:

```
make_bundle_from_anndata project.h5ad test_config_from_anndata.yaml test_bundle --write-premagetab
```

This will produce a subfolder in the output bundle called 'mage-tab' with a starting samples table, and in the case of droplet experiments a starting .cells.txt with the cell-wise information, which is distinct from sample-wise information for experiments of that type.

Pass these starting files to SCXA curators with all the experiment-level data you have. They will curate the metadata, in the process generating:

 - Accession the experiment (give it an identifier)
 - Create an SDRF (and cells.txt for droplet experiments)
 - Create an IDF (experiment description)

### 3b. Second run, with completed MAGE-TAB

When the curators are done, you can run the command again to produce the final bundle:

```
make_bundle_from_anndata atlas_anndata/data/E-MTAB-6077.project.h5ad test_config_from_anndata.yaml test_bundle --magetab-dir <path to mage-tab dir> --exp-name <new accession>
```

This bundle directory should now be ready to ingest into Single Cell Expression Atlas.


