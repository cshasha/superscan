print('.initializing')

import argparse
import sys, os
import pandas as pd
import numpy as np
import scanpy as sc
from joblib import load
from scipy.sparse import csr_matrix, csc_matrix
from scipy.stats import entropy
from tabulate import tabulate
import warnings; warnings.simplefilter('ignore')

parser = argparse.ArgumentParser()
parser.add_argument('--dataset', required=True, help='path to the dataset file')
parser.add_argument('--out', default='predictions', help='output filename')
opt = parser.parse_args()

features = pd.read_csv('features.csv',index_col=0).index
model = load('model.joblib')
    
print('..loading data')

if not os.path.exists(opt.dataset):
    sys.exit("The dataset file does not exist.")

if (opt.dataset[-4:] != '.csv') & (opt.dataset[-5:] != '.h5ad'):
    sys.exit('Error: Dataset format must be either csv or h5ad.')

    
if opt.dataset[-4:] == '.csv':
    data = pd.read_csv(opt.dataset, index_col=0)
    
    if 'CD4' in data.index:
        data = data.T

    overlap = list(set(data.columns) & set(features))
    
    if len(overlap) < 100:
        sys.exit('Error: Not enough feature overlap.')
        
    extra = list(np.setdiff1d(features,overlap))
    data = data[overlap]
    data[extra] = np.nan
    data = data[features]
        
if opt.dataset[-5:] == '.h5ad': 
    data = sc.read(opt.dataset)
    colnames = data.var.index
    rownames = data.obs.index
    
    overlap = list(set(colnames) & set(features))
    if len(overlap) < 100:
        sys.exit('Error: Not enough feature overlap.')
    
    if (type(data.X) == csc_matrix) | (type(data.X) == csr_matrix):
        data = data[:,overlap]
        data = pd.DataFrame.sparse.from_spmatrix(data.X)
    else:
        data = data[:,overlap]
        data = pd.DataFrame(data.X)
        
    data.columns = overlap
    data.index = rownames
    extra = list(np.setdiff1d(features,overlap))
    data[extra] = np.nan
    data = data[features]
    
print('...predicting') 

#add l2, l3 labels

Y_pred = model.predict(data)

print('....calculating entropy scores') 

probs = model.predict_proba(data)
ent = entropy(probs, base=2,axis=1)

cell_type = []
ent_score = []
counts = []
for x in pd.Series(Y_pred).unique():
    cell_type.append(x)
    counts.append(np.sum(Y_pred == x))
    ent_score.append(np.mean(ent[Y_pred == x]))

print(tabulate(set(zip(cell_type, counts, ent_score)), headers = ['cell type', 'number of cells', 'mean entropy']))
print('.....exporting results')
#add index
pd.DataFrame(data={'cell':data.index,'prediction':Y_pred, 'entropy score': ent}).to_csv(opt.out + '_superscan.csv',index=False)
print('......done')
