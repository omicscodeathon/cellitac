## Dataset

cellitac was developed and tested on the **10x Genomics PBMC Multiome
(unsorted, 10k)** reference dataset — a paired scATAC-seq + scRNA-seq assay
on cryopreserved peripheral blood mononuclear cells from a healthy human
donor, generated with the Chromium Next GEM Single Cell Multiome ATAC + Gene
Expression kit and sequenced on an Illumina NovaSeq 6000.

- **Source:** 10x Genomics (donor material from AllCells)
- **Alignment:** GRCh38 (hg38)
- **License:** Creative Commons Attribution 4.0 International (CC BY 4.0)
- **10x Genomics dataset page:**
  <https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-no-cell-sorting-10-k-1-standard-2-0-0>

### Files required by the pipeline

cellitac detects files by extension, so **file names do not matter**. From
the 10x page above, download the three files with these extensions into one
folder:

### Required Files

| File | Extension |
| --- | --- |
| [Filtered feature-barcode matrix](https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_unsorted_10k/pbmc_unsorted_10k_filtered_feature_bc_matrix.h5) | `.h5` |
| [ATAC fragments](https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_unsorted_10k/pbmc_unsorted_10k_atac_fragments.tsv.gz) | `.tsv.gz` |
| [Fragments tabix index](https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_unsorted_10k/pbmc_unsorted_10k_atac_fragments.tsv.gz.tbi) | `.tsv.gz.tbi` |

### Download via CLI

Alternatively, you can download these files directly from the command line:

```bash
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_unsorted_10k/pbmc_unsorted_10k_filtered_feature_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_unsorted_10k/pbmc_unsorted_10k_atac_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_unsorted_10k/pbmc_unsorted_10k_atac_fragments.tsv.gz.tbi

```

### Usage

Once the files are downloaded, run the following command to process them:

```bash
cellitac --input /path/to/downloaded/files --output /path/to/results

```



Additional datasets used for further testing are documented in the
[`accessions/`](../accessions/) directory.
