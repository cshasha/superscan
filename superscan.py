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
parser.add_argument('--out_prefix', default='predictions', help='output filename prefix')
opt = parser.parse_args()

features = pd.read_csv('features.csv',index_col=0)['0']
model_l1 = load('model_l1.joblib')
model_l2 = load('model_l2.joblib')
    
print('..loading data')

if not os.path.exists(opt.dataset):
    sys.exit("Error: The dataset file does not exist.")

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

Y_pred1 = model_l1.predict(data)
Y_pred2 = model_l2.predict(data)

print('....calculating entropy scores') 

probs1 = model_l1.predict_proba(data)
ent1 = entropy(probs1, base=2,axis=1)
ent1 /= np.log2(len(probs1[0]))
con1 = np.repeat('high', len(ent1))
con1[ent1 > 0.1] = 'med'
con1[ent1 > 0.3] = 'low'

probs2 = model_l2.predict_proba(data)
ent2 = entropy(probs2, base=2,axis=1)
ent2 /= np.log2(len(probs2[0]))
con2 = np.repeat('high', len(ent2))
con2[ent2 > 0.1] = 'med'
con2[ent2 > 0.3] = 'low'

cell_type = []
ent_score = []
counts = []
for x in pd.Series(Y_pred1).unique():
    cell_type.append(x)
    counts.append(np.sum(Y_pred1 == x))
    ent_score.append(np.mean(ent1[Y_pred1 == x]))

print('.....exporting results')

pd.DataFrame(data={'cell':data.index,'broad label':Y_pred1, 'normalized entropy score, broad': ent1, 'confidence, broad': con1,
                   'fine label':Y_pred2, 'normalized entropy score, fine': ent2,'confidence, fine': con2}).to_csv(opt.out_prefix + '_superscan.csv',index=False)
print('......done. Results saved to '+opt.out_prefix+'_superscan.csv')
print(tabulate(set(zip(cell_type, counts, ent_score)), headers = ['cell type', 'number of cells', 'mean entropy']))
