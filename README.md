# Single-cell Expression Atlas - compatible analysis bundles from annData files

In the Gene Expression team at the EBI we produce [Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/home), using a consistent pipeline to analyse data from the raw FASTQ files and produce the results visibile in the SCXA interface. Intermediate in this process is an 'analysis bundle' whereby the subset of results needed are gathered and formatted correctly for loading into our databases and indices. 

We sometimes encounter datasets without raw data, but where the user has an [annData](https://anndata.readthedocs.io/en/latest/) file containing expression matrices and derived results. The purpose of this repository is to provide a way of producing an analysis bundle directly from an annData object, allowing us (with some manual curation) to input these experiments. It will also be useful to simplify our own processes, since we can take the annData files now produced at the end of our analysis to produce a bundle in one step.

## Example commands  
