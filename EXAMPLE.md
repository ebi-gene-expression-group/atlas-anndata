# Example walkthrough for pre-analysed data in anndata format

This is an example walkthrough for processing externally analysed data for Single Cell Expression Atlas. For a more generic desrcription of tools and options as well as install instructions see the [README](README.md).

## Dataset

We'll we working with an annData file available from the GTEx portal [here](https://gtexportal.org/home/datasets) ([dataset link](https://storage.googleapis.com/gtex_analysis_v9/snrna_seq_data/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad)).

## Step 1: produce a starting config

```
make_starting_config_from_anndata --droplet GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad gtex_sc.yaml
```

This will produce a starting point of a config file that will eventually be used to make a final SCXA bundle. It has a number of 'FILL ME' tags indicating things that need your input.

## Step 2: first bundle run: write basic info

Without any editing of the config file we can output the basic info:

```
make_bundle_from_anndata --write-premagetab GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad gtex_sc.yaml gtex_sc_bundle
``` 

The output is like:

```
tree gtex_sc_bundle/
gtex_sc_bundle/
├── clusters_for_bundle.txt
├── mage-tab
│   ├── NONAME.precells.txt
│   └── NONAME.presdrf.txt
├── NONAME.cell_metadata.tsv
├── NONAME.project.h5ad
├── pca.tsv
├── reference
│   └── gene_annotation.txt
├── umap_tissue.tsv
├── umap.tsv
├── vae_mean_tissue.tsv
├── vae_mean.tsv
├── vae_samples.tsv
└── vae_var.tsv

2 directories, 13 files
```

## Step 3: Look at cellmetadata and refine configuration for curation

### Sample field

The sample field is encoded in the config generated above like:

```
  sample_field: FILL ME with a string
```

We need to identify if possible a field in the cell metadata that differentiates cells from different library preparations. If this field is not present this package will attempt to derive it from the cell IDs (in which case just unset the `FILL ME with a string` string in that field), which is what has happened by default in this case. Here are the first few lines of NONAME.cell_metadata.tsv, which was output above, and is just a straight dump of the cell metadata embedded in the object:

```
                                    id n_genes fpr         tissue prep
 CST01_TAGGCATGTAAATACG-skeletalmuscle    2658 0.1 skeletalmuscle  CST
 CST01_CCTTACGTCCGTCAAA-skeletalmuscle    2770 0.1 skeletalmuscle  CST
 CST01_CACACTCCATGGAATA-skeletalmuscle    2858 0.1 skeletalmuscle  CST
 CST01_CCACTACGTGCAACGA-skeletalmuscle    2899 0.1 skeletalmuscle  CST
 CST01_GCAGTTAAGACTGTAA-skeletalmuscle    2903 0.1 skeletalmuscle  CST
 CST01_CCTCTGAGTATCGCAT-skeletalmuscle    2484 0.1 skeletalmuscle  CST
 individual nGenes nUMIs PercentMito  PercentRibo Age_bin  Sex
          1   2902 11544 0.076230080 0.0511954280   51-60 Male
          1   3187  7727 0.004529572 0.0028471593   51-60 Male
          1   3344  7058 0.000000000 0.0009917824   51-60 Male
          1   3213  6487 0.000000000 0.0201942340   51-60 Male
          1   3251  6142 0.000000000 0.0244220120   51-60 Male
          1   2895  6119 0.000000000 0.0017976793   51-60 Male
                Sample.ID Participant.ID Sample.ID.short
 GTEX-1HSMQ-5011-SM-GKSJH     GTEX-1HSMQ GTEX-1HSMQ-5011q
 GTEX-1HSMQ-5011-SM-GKSJH     GTEX-1HSMQ GTEX-1HSMQ-5011
 GTEX-1HSMQ-5011-SM-GKSJH     GTEX-1HSMQ GTEX-1HSMQ-5011
 GTEX-1HSMQ-5011-SM-GKSJH     GTEX-1HSMQ GTEX-1HSMQ-5011
 GTEX-1HSMQ-5011-SM-GKSJH     GTEX-1HSMQ GTEX-1HSMQ-5011
 GTEX-1HSMQ-5011-SM-GKSJH     GTEX-1HSMQ GTEX-1HSMQ-5011
 RIN.score.from.PAXgene.tissue.Aliquot RIN.score.from.Frozen.tissue.Aliquot
                                   8.5                                  7.8
                                   8.5                                  7.8
                                   8.5                                  7.8
                                   8.5                                  7.8
                                   8.5                                  7.8
                                   8.5                                  7.8
 Autolysis.Score Sample.Ischemic.Time..mins. Tissue.Site.Detail scrublet
            None                         387  Muscle - Skeletal    False
            None                         387  Muscle - Skeletal    False
            None                         387  Muscle - Skeletal    False
            None                         387  Muscle - Skeletal    False
            None                         387  Muscle - Skeletal    False
            None                         387  Muscle - Skeletal    False
 scrublet_score          barcode batch n_counts tissue.individual.prep
      0.1840000 TAGGCATGTAAATACG    73     9745  skeletalmuscle_01_CST
      0.3513514 CCTTACGTCCGTCAAA    73     6230  skeletalmuscle_01_CST
      0.2968750 CACACTCCATGGAATA    73     5771  skeletalmuscle_01_CST
      0.1638142 CCACTACGTGCAACGA    73     5691  skeletalmuscle_01_CST
      0.1955307 GCAGTTAAGACTGTAA    73     5382  skeletalmuscle_01_CST
      0.1734694 CCTCTGAGTATCGCAT    73     5147  skeletalmuscle_01_CST
                   Broad.cell.type                 Granular.cell.type introns
 Myocyte (sk. muscle, cytoplasmic) Myocyte (slow-twitch, cytoplasmic)    9843
 Myocyte (sk. muscle, cytoplasmic) Myocyte (slow-twitch, cytoplasmic)   12430
                Myocyte (NMJ-rich)                 Myocyte (NMJ-rich)      NA
       Endothelial cell (vascular)      Endothelial cell (vascular) I   11675
       Endothelial cell (vascular)      Endothelial cell (vascular) I   10985
              Myocyte (sk. muscle)               Myocyte (sk. muscle)   14030
 junctions exons sense antisense intergenic
       467  6055 13989      2376     111339
       148  1703  7811      6470     102444
        NA    NA    NA        NA         NA
       189  1490  6928      6426      90757
       177  1343  6375      6130      82289
       124  1062  6034      9182     103031
                          batch.barcode exon_ratio intron_ratio junction_ratio
 skeletalmuscle_01_CST-TAGGCATGTAAATACG 0.36999694    0.6014665    0.028536511
 skeletalmuscle_01_CST-CCTTACGTCCGTCAAA 0.11924935    0.8703872    0.010363420
                                    nan         NA           NA             NA
 skeletalmuscle_01_CST-CCACTACGTGCAACGA 0.11157706    0.8742699    0.014153063
 skeletalmuscle_01_CST-GCAGTTAAGACTGTAA 0.10739704    0.8784486    0.014154338
 skeletalmuscle_01_CST-CCTCTGAGTATCGCAT 0.06979495    0.9220557    0.008149317
 log10_nUMIs leiden leiden_tissue Tissue.composition Cell.types.level.2
    4.062356     20             7             Muscle             Muscle
    3.888011     10             7             Muscle             Muscle
    3.848682     44            13             Muscle             Muscle
    3.812044      3             6       Endothelial    Endothelial cell
    3.788310      3             6       Endothelial    Endothelial cell
    3.786681     10             0             Muscle             Muscle
 Cell.types.level.3 Broad.cell.type.numbers
            Stromal                      33
            Stromal                      33
            Stromal                      29
            Stromal                       3
            Stromal                       3
            Stromal                      32
             Broad.cell.type..numbers.          Tissue
 33. Myocyte (sk. muscle, cytoplasmic) Skeletal muscle
 33. Myocyte (sk. muscle, cytoplasmic) Skeletal muscle
                29. Myocyte (NMJ-rich) Skeletal muscle
        3. Endothelial cell (vascular) Skeletal muscle
        3. Endothelial cell (vascular) Skeletal muscle
              32. Myocyte (sk. muscle) Skeletal muscle
                       channel sample
 skeletalmuscle_CST_GTEX-1HSMQ  CST01
 skeletalmuscle_CST_GTEX-1HSMQ  CST01
 skeletalmuscle_CST_GTEX-1HSMQ  CST01
 skeletalmuscle_CST_GTEX-1HSMQ  CST01
 skeletalmuscle_CST_GTEX-1HSMQ  CST01
 skeletalmuscle_CST_GTEX-1HSMQ  CST01
```

** Note: the 'sample' column on the end has been guessed from the cell IDs and was not a formal part of the cell metadata **

Referring to [the paper](https://www.biorxiv.org/content/10.1101/2021.07.19.452954v1.full), the guessed sample derived from the cell IDs appears to a combination of single-nucleus isolation protocol and individual. I'm not sure why the object only has 4 individuals, since that is not what the paper described. In any case, we can see that here the cell IDs also contained the tissue, so it seems sensible to use `tissue individual prep` as our sample field, so we edit the yaml accordingly:

```
  sample_field: tissue-individual-prep
```

### Flag curated fields

Cell meta data from annData objects is a mixtrue of any input sample metadata provided by the author, plus annotations added over the course of analysis. The latter may not be appropriate for inclusion in the metadata in SCXA. You can edit the `cell_metadata` section of gtex_sc.yaml:

```
...
cell_meta:
  entries:
  - default: false
    kind: analysis
    markers: false
    parameters: {}
    slot: n_genes
  - default: false
    kind: curation
    markers: false
    parameters: {}
    slot: fpr
...
```

... and change the `kind` to `analysis` for any field that should not be included in the `precells` and `presdrf`. 


## Step 4: second bundle run: re-write bundle starting metadata for curation

With the updated configuration in hand we can generate a bundle with pre-MAGE-TAB data for curation:

```
make_bundle_from_anndata --write-premagetab GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad gtex_sc.yaml gtex_sc_bundle
``` 

This will produce the pre-MAGE-TAB files at `/mage-tab/NONAME.presdrf.txt` in the bundle.

## Step 5: Do curation:

### Regular curation efforts

The pre-MAGE-TAB can now be used to start curation by the curation team.

### Gather other missing info needed before final bundle creation

Unlike our standard submission pathways, for pre-analysed data we need additional information before the data are ingested for SCXA, which must currently be provided via the configuration YAML. The completed config from the following steps should be added to the `scxa-metadata` alongside the MAGE-TAB files.

#### Identify correct gene ID and gene name fields

For SCXA we need the gene symbol and ID fields. **We only work with datasets keyed by Ensembl gene identifier**.

In the configuration step the package will try to guess the gene identifier and symbol fields. In this case we got:

```
gene_meta:
  id_field: index
  name_field: gene_name
```

** Note: 'index' is just the row names of the gene annotation data in the anndata object **

Have a look at the `reference/gene_annotation.txt` file from the bundle dir:

```
gene_id	gene_ids	Chromosome	Source	Start	End	Strand	gene_name	gene_source	gene_biotype	gene_length	gene_coding_length	Approved symbol	Approved nameStatus	Previous symbols	Alias symbols	gene_include
FO538757.2	ENSG00000279457	1	ensembl	184922	200322	-	FO538757.2	ensembl	protein_coding	15400	1982	WASH9P	WAS protein family homolog 9, pseudogene	Approved	nan	nan	True
SAMD11	ENSG00000187634	1	ensembl_havana	924879	944581	+	SAMD11	ensembl_havana	protein_coding	19702	3214	SAMD11	sterile alpha motif domain containing 11	Approved	nan	MGC45873	True
NOC2L	ENSG00000188976	1	ensembl_havana	944203	959309	-	NOC2L	ensembl_havana	protein_coding	15106	5539	NOC2L	NOC2 like nucleolar associated transcriptional repressor	Approved	nan	DKFZP564C186, NET7, NET15, NIR, PPP1R112	True
KLHL17	ENSG00000187961	1	ensembl_havana	960586	965715	+	KLHL17	ensembl_havana	protein_coding	5129	3395	KLHL17	kelch like family member 17	Approved	nan	nan	True
```

We can see that the guess has not worked here- the gene metadata clearly has gene symbols rather than IDs in its row names. Not great practice. In any case we need to edit the config to use the `gene_ids` field:

```
gene_meta:
  id_field: gene_ids
  name_field: gene_name
```

#### Information on reference and software versions etc

```
analysis_versions:
- analysis: reference
  asset: FILL ME with a string
  citation: FILL ME with a string
  kind: file
  version: FILL ME with a string
- analysis: filtering and trimming
  asset: FILL ME with a string
  citation: FILL ME with a string
  kind: software
  version: FILL ME with a string
- analysis: mapping
  asset: FILL ME with a string
  citation: FILL ME with a string
  kind: software
  version: FILL ME with a string
- analysis: clustering
  asset: FILL ME with a string
  citation: FILL ME with a string
  kind: software
  version: FILL ME with a string
```

As many of these as possible should be filled, any incomplete entries should be removed. For example we might set things like:

```
analysis_versions:
- analysis: reference
  asset: Ensembl
  kind: file
  version: '84'
- analysis: mapping
  asset: CellRanger
  kind: software
  version: various
```

#### Matrices

This is more important, since it defines the extent to which expression data from a submitted experiment can be used in the SCXA database. Essentially we need a normalised count matrix. All matrices will also be exported for providsion on our FTP. This config section looks like:

```
matrices:
  entries:
  - cell_filtered: true
    gene_filtered: false
    log_transformed: false
    measure: FILL ME with a string
    name: FILL ME with a string
    normalised: false
    parameters: {}
    scaled: true
    slot: raw.X
  - cell_filtered: true
    gene_filtered: false
    log_transformed: false
    measure: FILL ME with a string
    name: counts
    normalised: false
    parameters: {}
    scaled: true
    slot: counts
  - cell_filtered: true
    gene_filtered: false
    log_transformed: false
    measure: FILL ME with a string
    name: FILL ME with a string
    normalised: false
    parameters: {}
    scaled: true
    slot: X
  load_to_scxa_db: FILL ME with a string
```

The status flags for each matrix MUST be set correctly, and a name supplied (which will be used for the exported matrix file). Where information cannot be discovered, the whole entry for a matrix should be removed.


## Step 6: third bundle run: write final bundle from completed YAML file

```
make_bundle_from_anndata --exp-name <exp name> GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad gtex_sc.yaml --conda-prefix <conda prefix> gtex_sc_bundle
```

This will pull the curated metadata from the `scxa-metadata` repo, condense the SDRF (adding ontology terms) and re-generated cell-wise annotations that will be used to enrich the content of the annData file. This workflow uses conda to pull the relevant software dependencies.


 - `exp name` is the E- accession generated for the experiment in curation
 - `conda prefix` is location conda environments will be stored. 
