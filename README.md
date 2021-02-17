# Superscan: Supervised Single-Cell Annotation

## Overview

Superscan (Supervised Single-Cell Annotation) is a Python package for annotation of human scRNA-seq data. Superscan is a supervised model that was trained on a collection of publicly available datasets that were manually labeled using corresponding CITE-seq data.

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

Superscan is run directly from the command line. It can read data in either anndata (.h5ad) or csv format (where the csv contains a counts matrix, with genes as columns). Gene names must be gene symbols (not ENSEMBL IDs), and raw counts (not normalized) should be provided.

```
usage: python superscan.py [-h] --dataset DATASET [--out OUT]

required arguments:
  --dataset DATASET  path to the dataset file

optional arguments:
  -h, --help                show this help message and exit
  --out_prefix OUT_PREFIX   output filename prefix
  ```
  
  The only required argument is the input dataset. The output filename can be specified, and will default to "predictions". In addition to saving the results in this file, a summary of the predictions will be shown, displaying the number of cells classified to each cell type, as well as the mean normalized entropy score (an indication of the confidence of the prediction) per cell type. Normalized entropy scores less than 0.1 indicate high confidence in the classification, normalized entropy scores between 0.1 and 0.3 indicate medium confidence, and normalized entropy scores greater than 1 indicate low confidence.
  
  Example datasets can be found in the `test_datasets` folder. An example command, using a dataset found at `test_datasets/test_ds.h5ad`, would be:
  
  ```
  python superscan.py --dataset test_datasets/test_ds.h5ad --out_prefix test
  ```
  The resulting predictions as well as corresponding entropy scores will be saved to `test_superscan.csv`.
