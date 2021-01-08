# Superscan: Supervised Single-Cell Annotation

## Overview

Superscan (Supervised Single-Cell Annotation) is a Python package for annotation of huamn scRNA-seq data. Superscan is a supervised model that was trained on a collection of publicly available datasets that were manually labeled using corresponding CITE-seq data.

## Installation

Superscan requires Python 3. The required files can be downloaded by cloning this repository into the desired directory:

```
git clone https://github.com/cshasha/superscan
```

The package dependencies can then be installed using a virtual environment:

```
python3 -m venv superscan_env
source superscan_env/bin/activate
pip install -r requirements.txt
```

## Cell Annotation

Superscan is run directly from the command line. It can read data in either anndata (.h5ad) or csv format. Gene names must be gene symbols (not ENSEMBL IDs). 

```
usage: python superscan.py [-h] --dataset DATASET [--out OUT]

required arguments:
  --dataset DATASET  path to the dataset file

optional arguments:
  -h, --help         show this help message and exit
  --out OUT          output filename
  ```
  
  The only required argument is the input dataset. The output filename can be specified, and will default to "predictions". In addition to saving the results in this file, a summary of the predictions will be shown, displaying the number of cells classified to each cell type, as well as the mean entropy score (an indication of the confidence of the prediction) per cell type. Entropy scores less than 0.5 indicate high confidence in the classification, entropy scores between 0.5 and 1 indicate medium confidence, and entropy scores greater than 1 indicate low confidence.
  
  An example command, using a dataset named test_ds.h5ad, would be:
  
  ```
  python superscan.py --dataset test_datasets/test_ds.h5ad --out test
  ```
  The resulting predictions as well as corresponding entropy scores will be saved to test_superscan.csv.
