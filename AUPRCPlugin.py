#!/usr/bin/env python
# coding: utf-8

# ## AUPR for evaluating taxonomic profiles in the benchmarking work 
# * Dec. 22, 2020
# 
# This notebook is intented to illustrate how to calculate the AUPRC based on an expected abundance profile and observed one. 

# In[4]:


import pandas as pd
import numpy as np
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import auc
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
import matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')


# ### precision and recall in general
# * Precision: the ratio tp / (tp + fp) where tp is the number of true positives and fp the number of false positives. The precision is intuitively the ability of the classifier not to label as positive a sample that is negative.
# 
# * Recall: the ratio tp / (tp + fn) where tp is the number of true positives and fn the number of false negatives. The recall is intuitively the ability of the classifier to find all the positive samples.
# 

# ### precision and recall in the context of evaluating the performance of taxonomic profiles
# 
# * Precision: the proportion of true classified taxa over all classified taxa. 
# * Recall: the proportion of correctly classified abundances over all truth abundances
# 

# ### the updated function for calculating precision recall curve
# * we modified the recall=1.0 with the higest observed recall in the resulting recall vector,
# * It is also supported by the CELL paper, where they "set the precision to 0 for the highest observed recall onward" (https://www.cell.com/cms/10.1016/j.cell.2019.07.010/attachment/3a60e4f3-3b48-4909-8c4c-c3da3fd46267/mmc1.pdf)

# In[14]:


def to_labels(abd, threshold):
    return (abd>threshold).astype('int')

def updated_prc(expected_abd, observed_abd, thresholds=[0, 0.1, 0.001]):
    bi_expected_abd=np.array(expected_abd!=0, dtype=int)
    precision, recall, _ = precision_recall_curve(bi_expected_abd, observed_abd)
    precision[0]=max(precision[1:])
    recall[0]=max(recall[1:])
    aupr=auc(recall, precision)
    f1_scores = [f1_score(bi_expected_abd, to_labels(observed_abd, i)) for i in thresholds]
    precision_scores = [precision_score(bi_expected_abd, to_labels(observed_abd, i)) for i in thresholds]
    recall_scores = [recall_score(bi_expected_abd, to_labels(observed_abd, i)) for i in thresholds]
    return np.append(np.concatenate((precision_scores, recall_scores, f1_scores)), aupr)


# Example:

# In[16]:
class AUPRCPlugin:
   def input(self, infile):

    expected_abd=np.array([0, 0.4, 0.3, 0, 0.2, 0.1, 0])
    observed_abd=np.array([0.1, 0, 0.2, 0, 0.7, 0, 0])
    updated_prc(expected_abd, observed_abd)


    expected_abd=np.array([0, 0.4, 0.3, 0.2, 0.1, 0])
    observed_abd=np.array([0.1, 0, 0.2, 0.7, 0, 0])
    updated_prc(expected_abd, observed_abd)

    file=infile#"D_M_K_m_Raw_AbdTable_1225.txt"

    self.data=pd.read_csv(infile, sep='\t', index_col=0)

   def run(self):
    self.data[self.data['Environment']=='Building1']

    ids = self.data.Environment.unique().tolist()
    n_ids = len(ids)

    tmp=self.data[self.data['Environment']=='Building1']
    for i in range(2, 6):
        print(updated_prc(tmp.iloc[0, 3:], tmp.iloc[i, 3:]))
        print(updated_prc(tmp.iloc[1, 3:], tmp.iloc[i, 3:]))


    thresholds = np.concatenate((np.arange(0.0, 0.0001, 0.00001), 
                             np.arange(0.0001, 0.001, 0.0001), 
                             np.arange(0.001, 0.011, 0.001))).tolist()

    thresholds_str = ['{:.5f}'.format(x) for x in thresholds]
    thresholds_num =[float(x) for x in thresholds_str]


    metrics = ['precision', 'recall', 'f1']
    metrics_rep = [item for item in metrics for i in range(len(thresholds))]
    thresholds3 = thresholds_str*len(thresholds)
    concat_func = lambda x,y: x + "_" + str(y)

    name_all = list(map(concat_func, metrics_rep,  thresholds3))
    name_all.append('aupr')
    tax_name_all = list(map(concat_func, ["tax"]*len(name_all), name_all))
    seq_name_all = list(map(concat_func, ["seq"]*len(name_all), name_all))
    colnames = tax_name_all + seq_name_all


    ncols=(len(thresholds_num)*3+1)*2
    out = []
    rownames = []
    for i in ids: #xrange(len(ids)):
        tmp=self.data.groupby('Environment').get_group(i)
        rownames.extend(tmp.index.values[2:].tolist())
        #sample_arr = np.empty([4, ncols])
        for j in range(2, 6):
            tax = updated_prc(tmp.iloc[0, 2:], tmp.iloc[j, 2:], thresholds_num)
            seq = updated_prc(tmp.iloc[1, 2:], tmp.iloc[j, 2:], thresholds_num)
            out.append(np.concatenate((tax, seq)))
    np_out=np.array(out)

    self.out=pd.DataFrame(data=np_out, columns=colnames, index=rownames)
    self.out

   def output(self, outputfile):
      self.out.to_csv(outputfile, sep='\t')


