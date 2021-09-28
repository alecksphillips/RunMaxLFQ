# RunMaxLFQ
Simple script for applying MaxLFQ algorithm post-hoc to MaxQuant output with
option to distribute batches of protein groups across cores/cluster

Uses implementations of MaxLFQ from [iq](https://cran.r-project.org/web/packages/iq/index.html) and [diann](https://github.com/vdemichev/diann-rpackage)

## Requirements

### CRAN Packages
* data.table
* bit64
* fst
* progress
* pbapply
* parallel
* iq

### Github R Packages
* [diann](https://github.com/vdemichev/diann-rpackage)

## Usage

```R
Rscript --vanilla RunMaxLFQ.R \
  <path to evidence/peptides.txt> \
  <path to proteinGroups.txt, default: ./proteinGroups.txt> \
  <path for output, default: ./output.txt> \
  <serial/parallel, default: serial> \
  <nBatches default: 1> <nWorkers, default: 1> \
  <useDIANN, default: FALSE>
```