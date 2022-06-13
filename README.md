[![Anaconda-Server Badge](https://anaconda.org/ebi-gene-expression-group/atlas-anndata/badges/installer/conda.svg)](https://anaconda.org/ebi-gene-expression-group/atlas-anndata) [![PyPI version fury.io](https://badge.fury.io/py/atlas-anndata.svg)](https://pypi.python.org/pypi/atlas-anndata/)
[![Build Status](https://api.travis-ci.com/ebi-gene-expression-group/atlas-anndata.svg?branch=develop)](https://travis-ci.org/ebi-gene-expression-group/atlas-anndata)


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

The steps required to produce a bundle are:

 1. Initialise the bundle (`init` step). Process the annData file to determine what information is available and tore a summary of this information in a YAML-format configuration file (`anndata-config.yaml` in the bundle directory), and us that to Generate bundle files (1st time) based on the starting config files in 1. 
 2. Examine the cell metadata files and use that information to refine configuration related to cell metadata. This includes flagging fields that should be included in the pre-MAGE-TAB files, and for droplet experiments, finding the field that separates cells from different libraries. Also check the gene metadata at this stage, to check that Ensembl gene identifiers are available.
 3. Initialise MAGE-TAB (`init_magetab` step). Generate bundle files (2nd time) based on configuration refined in 2.. This will include pre-MAGE-TAB files suitable as a basis for curation.
 4. Undertake curation based on pre-magetab files generated at 3.
 5. Inject curated metadata into the annData object and configuration (`inject_magetab` step). Condense the new SDRF file (including zoomification) and inject the new metadata back into the annData object where it can be used in curation. Output a new configuration that includes any new fields added in curation.
 6. Make final edits to the config YAML. Complete anything marked with FILL ME, set the `load_to_scxa_db` matrix and flag any e.g. cell type fields for marker generation.
 7. Generate the final bundle suitable for loading into SCXA (`final` step).

Once finalised, the config YAML should be added to the `scxa_metadata` repo alongside the MAGE-TAB files.

## Role of curators and bioinformaticians

The process is progressive, with information on an experiment gradulally being accumulated for potential partial analysis and loading. Steps 1-5 could be undertaken by curators, with the bundle and configuration then being handed off to bionformaticians. However additional author-provided information (e.g. on matrix processing) is likely to be needed at step 6, and it may be that that information could be gathered from authors alongside curation in step 4.

## Protocol: producing an analysis bundle from arbitrary annData files

### 0. Install this package

This repository contains a python packages which should be used to facilitate production of a bundle. 

Probably the simplest install method is with Conda (or Mamba). e.g. to create a new environment with the package:

```
conda create -n atlas-anndata atlas-anndata
```

OR install it from PyPi like:

```
pip install atlas-anndata
```

The `make_bundle_from_anndata` script is the central utility over which a bundle is created over the steps above. The usage is:

```
Usage: make_bundle_from_anndata [OPTIONS] EXP_NAME
                                [init|init_magetab|inject_magetab|final]

  Build a bundle directory compatible with Single Cell Expression Atlas
  (SCXA) build proceseses

  exp_name       - Specify an Expression Atlas identifier for this experiment.
  step           - Specify the bundle creation step. One of:
                      * 'init': start a bundle based on anndata file content
                      * 'init_magetab': create starting MAGE-TAB files as a
                        basis for curation, having checked and refined
                        configuration at the 'init' stage.
                      * 'inject_magetab': With curation done, read metadata
                        from the scxa-metadata repo, and modify bundle
                        configuration and annData object accordingly.
                      * 'final': Having made any refinements to the field
                        configuration modified by 'inject_magetab', produce the
                        final bundle.

Options:
  --anndata-file PATH            For the 'init' stage, specify a path to a
                                 file of the  annData hdf5 specification, with
                                 all necessaryinformation for SCXA.

  --bundle_dir PATH              A directory under which bundle directories
                                 should be created. Defaults to the current
                                 working directory.

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

  --max-rank-for-stats INTEGER   For how many top marker genes should stats
                                 (mean, median expression) be output?

  --matrix-for-markers TEXT      Where cell groups in the configuration file
                                 have been flagged with markers, which matrix
                                 should be used? Can be X, or an entry in
                                 .layers(). The matrix must be appropriate for
                                 Scanpy's tl.rank_genes_groups() method,
                                 usually meaning filtered, normalised and log
                                 transformed, but without additional scaling.
                                 [required]

  --conda-prefix PATH            Specify a Conda directory to be used for
                                 environments when running Snakemake
                                 workflows.

  --scxa-metadata-branch TEXT    When searching the SCXA metadata repository
                                 for curation for this experiment, which
                                 branch should we use?  [required]

  --sanitize-columns             When adding data from curation into the
                                 anndata object, should we remove the Comment,
                                 Characteristic etc?

  --exp-name TEXT                Specify an Expression Atlas identifier that
                                 will be used for this experiment. If not set,
                                 a placeholder value E-EXP-1234 will be used
                                 and can be edited in the bundle later.
                                 [required]

  --scxa-db-scale INTEGER        To what overall scale should cell counts be
                                 multiplied for the SCXA DB? A multiplier will
                                 be calculated from this value and the median
                                 cell-wise sum in the given matrix.

  --help                         Show this message and exit.
```

(Note that the `atlas-style` flag is probably only useful for annData files produced by the Experession Atlas team, and relies on a number of assumptions about the content of the file in order to infer some additional information.)

### 1. Bundle initialisation: Produce a YAML format annData description file and starting bundle content

To produce a valid bundle from an anndata file, we need to describe that file, outlining which of the cell/ gene metadata columns, matrices,dimension reductions etc should be included. This is done via a YAML-format config file (see [example](atlas_anndata/data/bundles/E-MTAB-6077/anndata-config.yaml)). Then produce a starting bundle based on that config. 

Both of these things are accomplished by 'make_starting_config_from_anndata` with the 'init' step:

```
make_bundle_from_anndata --anndata-file atlas_anndata/data/bundles/E-MTAB-6077/E-MTAB-6077.project.h5ad E-MTAB-6077 init
``` 

(note: we would supply `--droplet` at this stage for a droplet experiment).

This will create a starting bundle, by default in the current directory (see `--bundle-dir`), including a starting version of the bundle configuration and a copy of the annData file which we will eventually customise and distribute:

```
> tree E-MTAB-6077/
E-MTAB-6077/
├── anndata-config.yaml
├── E-MTAB-6077.cell_metadata.tsv
├── E-MTAB-6077.project.h5ad
├── MANIFEST
└── reference
    └── gene_annotation.txt

1 directory, 5 files
```

The config is likely to wrong in a number of ways, but its just a starting point. 

### 2: refine configuration for curation

With the configuration file and unmoderated content available to us we can make some sensible decisions about some of those settings in the YAML file.

#### Identify correct gene ID and gene name fields

For SCXA we need the gene symbol and ID fields. The configuration YAML might have populated these if the default field names are present, but you may well get:

```
gene_meta:
  id_field: FILL ME with a string
  name_field: FILL ME with a string
```

You need to look at the `reference/gene_annotation.txt` file from the bundle directory and set these fields. `id_field` **must** be a field containing Ensembl gene IDs. If these are not available we cannot work with a dataset. `name_field` is a field containing gene symbols.

### Sample field (droplet experiments only)

The sample field is encoded in the config generated above like:

```
  sample_field: FILL ME with a string
```

The value of this configuration field must be a field name from NONAME.cell_metadata.tsv in the bundle directory corresponding to a field that separates cells from different libraries.The `sample` field in this file is usually derived from the cell identifiers, and should only be used in the absence of more concrete information.

### 3: Initialise MAGE-TAB

Starting MAGE-TAB content can now be generated with the `init_magetab` step of `make_bundle_from_anndata`:

```
make_bundle_from_anndata E-MTAB-6077 init_magetab
```

The experiment name ('E-MTAB-6077' here) and step ('init_magetab') must be supplied. The annData and configuration will be read from the bundle directory, assumed by default to be in the current working directory (see above).

### 4: Do curation

The pre-MAGE-TAB can now be used to start curation by the curation team. Assuming contact with authors is occurring at this point, it would save time if the information needed at step 6 could also be gathered here, such as the processing status of component matrices. Otherwise the appropriate contact should be passed to bioinformaticians, who will need this info prior to producing the final bundle.

### 5. Inject curated metadata 

MAGE-TAB format metadata should now be available in the `scxa-metadata` repository, from where we can retrieve it to be condensed, zoomified, and added to the annData object. This can all be done with another run of `make_bundle_from_anndata`:

```
make_bundle_from_anndata E-MTAB-6077 inject_magetab
```

This will pull the curated metadata from the `scxa-metadata` repo, condense the SDRF (adding ontology terms) and re-generate cell-wise annotations that will be used to enrich the content of the annData file, with any new fields from curation added to `.obs`. It will also add configuration for those fields to the YAML.

Not that `--conda-prefix` can be spedified here, and is a location conda environments will be stored for the workflow that does SDRF condensation etc.

### 6. Finalise configuration

We're almost ready to creat our final bundle, but we must first finalise all the configuration in `anndata-config.yaml` for the bundle.

#### Flag curated fields

Cell meta data from annData objects is a mixture of any input sample metadata provided by the author, plus annotations added over the course of analysis. The latter may not be appropriate for inclusion in the metadata in SCXA. Check the fields described in `cell_meta`, especially their kind ('curation', 'clustering', 'analysis'). Curated fields are those present before analysis, biological and technical info for cells and samples. Clustering is used to indicate the results of unsupervised cell clustering stored in .obs. Analysis is everthing else, comprising all other fields added to .obs during analysis.

Most importantly: 

 - Flip `curation` to `analysis` for any field entry which should not ultimately form part of the SCXA experiment MAGE-TAB format metadata.
 - Ensure any fields corresponding to unsupervised clusterings are flagged correctly.

#### Gather other missing info needed before final bundle creation

Unlike our standard submission pathways, for pre-analysed data we need additional information before the data are ingested for SCXA, which must currently be provided via the configuration YAML. The completed config from the following steps should be added to the `scxa-metadata` alongside the MAGE-TAB files.

 - Under analyses please describe the analysis that was done. At a minimum you should describe the reference used (see [the example](atlas_anndata/example_config.yaml)) and the mapping tool used.
 - Under `matrices` check that you want all these matrices to be considered. You can remove any matrix that's not useful, and you should check the processing flags / matrices for each one. 
 - Under load_to_scxa_db please state the matrix that should be used by Atlas in expression-based displays. This should be filtered and normalised but not scaled or transformed. If no matrix in the object matches these criteria please remove this part of the config and Atlas will not show displays for this experiment based on expression values.
 - Check the dimension reductions described, again paying attention to 'kind'. 

For all sections, check [the example](atlas_anndata/data/bundles/E-MTAB-6077/anndata-config.yaml) for an idea of how things should look. For example under `matrices` there will be an entry pertaining the content of .X. You shold add a name (e.g 'scaled'), and check that the status of all the boolean flags is correct. This may well require some queries to the authors.

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

### Mark cell metadata fields for marker detection

Any categorical field can be used for marker detection (where appropriate matrices are available). This involves flipping the 'marker' status on a field annotation, e.g. changing:

```
  - default: false
    kind: clustering
    markers: true
    parameters: {}
    slot: leiden
```

... to 

```
  - default: false
    kind: clustering
    markers: false
    parameters: {}
    slot: leiden
```

#### Validate config YAML

Having edited the config YAML, you should validate it against a schema we provide and the annData file itself. We can use this mechanism to ensure that inputs match the expectations of Single Cell Expression Atlas. 

```
validate_anndata_with_config E-MTAB-6077
```

This will check that the config in the bundle directory for this experiment matches with the modified annData object there.

Bundling steps will also run this automatically before proceeding, but running it yourself will flag any issues early. If the validation flags any issues, resolve them.


### 7: Final bundle run

With the configuration finalised, it's just a case of running the `final` step for `make_bundle_from_anndata`:

```
make_bundle_from_anndata E-MTAB-6077 final
```

All bundle content including matrices, dimension reductionss etc will be written, and this should form the final bundle that can be read by SCXA loading processes.

