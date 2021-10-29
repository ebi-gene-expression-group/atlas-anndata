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

## Command help

Here is the description of options for the main bundle generation script:

```
> bundle_from_anndata.py --help
Usage: bundle_from_anndata.py [OPTIONS] ANNDATA_FILE BUNDLE_DIR

  Build a bundle directory compatible with Single Cell Expression Atlas (SCXA)
  build proceseses

  anndata_file - A file of the annData hdf5 specification, with all necessary information for SCXA.
  bundle_dir   - A directory in which to create the bundle.

Options:
  --exp-name TEXT                 Specify an Expression Atlas identifier that
                                  will be used for this experiment. If not
                                  set, a placeholder value E-EXP-1234 will be
                                  used and can be edited in the bundle later.
  --droplet                       Is this a droplet experiment?
  --exp-desc PATH                 Provide an experiment description file. If
                                  unspecified, the script will check for a
                                  slot called "experiment" in the .uns slot of
                                  the annData object, and use that to create a
                                  starting version of an IDF file.
  --write-cellmeta / --no-write-cellmeta
                                  Write a table of cell-wise metadata from
                                  .obs.
  --write-genemeta / --no-write-genemeta
                                  Write a table of gene-wise metadata from
                                  .var.
  --nonmeta-obs-patterns TEXT[,TEXT...]
                                  A comma-separated list of patterns to be
                                  used to exlude columns when writing cell
                                  metadata from .obs. Defaults to
                                  louvain,n_genes,n_counts,pct_,total_counts
  --nonmeta-var-patterns TEXT[,TEXT...]
                                  A comma-separated list of patterns to be
                                  used to exlude columns when writing gene
                                  metadata from .obs. Defaults to mean,counts,
                                  n_cells,highly_variable,dispersion
  --raw-matrix-slot TEXT          Slot where the most raw and least filtered
                                  expression data are stored. This is usually
                                  in raw.X, other values will be interpreted
                                  as layers (in .raw if prefixed with "raw".
                                  annData requires that .raw.X and .X match in
                                  .obs, so even when stored in .raw this will
                                  have to be the matrix afer cell filtering,
                                  but should ideally not be gene-filtered.
  --filtered-matrix-slot TEXT     Layer name or "X", specifying storage
                                  location for gene-filtered matrix before
                                  transformations such as log transform,
                                  scaling etc.
  --normalised-matrix-slot TEXT   Layer name or "X", specifiying storage
                                  location for normalised matrix.
  --final-matrix-slot TEXT        Layer name or "X", specifiying storage
                                  location for the final matrix after
                                  processing. Will usually be in .X.
  --write-obsms / --no-write-obsms
                                  Write dimension reductions from .obsm to the
                                  bundle?
  --obsms TEXT[,TEXT...]          A comma-separated list of obsm slots to
                                  write. Default is to write them all.
  --write-clusters / --no--write--clusters
                                  Write cluster data to the bundle?
  --clusters TEXT[,TEXT...]       A comma-separated list of .obs slots
                                  corresponding to unsupervised clusterings.
  --clusters-field-pattern TEXT   If --clusters not supplied, a string to be
                                  used to select columns from .obs
                                  representing unsupervised clusterings.
  --default-clustering TEXT       Where --write-clusters is set, which
                                  clustering should be set as the default? If
                                  not set, the middle (or first middle)
                                  clustering will be selected, or if --atlas-
                                  style is set, this will be the clustering
                                  corresponding to a resolution of 1.
  --write-markers / --no--write--markers
                                  Write marker data to the bundle?
  --marker-clusterings TEXT[,TEXT...]
                                  A comma-separated list of clusterings for
                                  which to write markers. marker results are
                                  expected to be stored in .uns under keys
                                  like "markers_{clustering}. Defaults to all
                                  selected clusterings.
  --metadata-marker-fields TEXT[,TEXT...]
                                  A comma-separated list of .obs slots
                                  corresponding to metadata variables for
                                  which markers have been derived.
  --write-marker-stats / --no--write--marker-stats
                                  Write marker summary statistics (mean,
                                  median expression) to the bundle?
  --marker-stats-layers TEXT[,TEXT...]
                                  A comma-separated list of layers from which
                                  expression values should be summarised in
                                  marker statistics.
  --max-rank-for-stats INTEGER    For how many top marker genes should stats
                                  (mean, median expression) be output?
  --atlas-style TEXT              Assume the tight conventions from SCXA, e.g.
                                  on .obsm slot naming
  --gene-name-field TEXT          Field in .var where gene name (symbol) is
                                  stored.
  --write-anndata / --no--write--anndata
                                  Write the annData file itself to the bundle?
  --help                          Show this message and exit.
```

## Example commands  
