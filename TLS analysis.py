import pandas as pd
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster, ward
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
from statsmodels.stats.multitest import multipletests
#tls=pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\TLS marker gene.csv")
#gse113=pd.read_excel(r"C:\Users\Heather P\Desktop\github\T1\GSE119336_RNA.xlsx")
#gse107=pd.read_excel(r"C:\Users\Heather P\Desktop\github\T1\GSE107943_RPKM.xlsx")
#tcga=pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\result_df_60483gene.csv")
#tls=tls.iloc[:7,:]
#tcga = tcga.loc[:, tcga.columns.isin(tls['ensembl'])]
#gse113 = gse113.loc[gse113['Gene'].isin(tls['gene']),:]
#gse107 = gse107.loc[gse107['Ensenble'].isin(tls['ensembl']),:]#correct
#gse107=gse107.iloc[:,2:]
#############################################
tcga=pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\result_df_TLS.csv")
gse107=pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\GSE107943_TLS.csv")
gse119=pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\GSE119336_TLS.csv")
def log2w0(value):
    if value != 0:
        return np.log2(value)
    else:
        return value
#tcga.iloc[:, 1:8] = tcga.iloc[:, 1:8].applymap(log2w0)
#tcga=tcga.iloc[:,1:]
#sorting_col= 'Ward_D'
#tcga= tcga.sort_values(by=sorting_col)
#tumor_group=tcga.iloc[:,7]
#tcga=tcga.T

#x=gse107.iloc[:7,1:]
#for column in x.columns:
#    x[column] = pd.to_numeric(x[column])
#x=x.applymap(log2w0)
#gse107.iloc[:7, 1:] = x
#gse107.index=gse107.iloc[:,0]
#gse107=gse107.iloc[:,1:]
#gse107=gse107.T
#sorting_col= 'group'
#gse107= gse107.sort_values(by=sorting_col)
#gse107=gse107.T
#tumor_group=gse107.iloc[7,:]

x=gse119.iloc[:7,1:]
for column in x.columns:
    x[column] = pd.to_numeric(x[column])
x=x.applymap(log2w0)
gse119.iloc[:7, 1:] = x
gse119.index=gse119.iloc[:,0]
gse119=gse119.iloc[:,1:]
gse119=gse119.T
sorting_col= 'group'
gse119= gse119.sort_values(by=sorting_col)
#gse119=gse119.T
#tumor_group=gse119.iloc[7,:]
#hot_group, cold_group = 1, 2
#cluster_colors = {
#    hot_group: "#FDFEBA",
#    cold_group: "#00BBCC"}
#col_colors = np.array([cluster_colors[label] for label in tumor_group])
#sns.clustermap(gse119, col_cluster=False,row_cluster=False, z_score=0,cmap="bwr", col_colors=col_colors)
#plt.show()
#################################Compute p values
#######computation
def pv(df1):
    res = pd.DataFrame()
    hot = gse119[gse119['group'] == 1]
    cold = gse119[gse119['group'] == 2]
    for i in df1.columns[:-1]:
        u_statistic, p_value = stats.mannwhitneyu(hot[i], cold[i], alternative='two-sided')
        result_df = pd.DataFrame({'MannWhitneyU_Statistic': [u_statistic], 'pvalue': [p_value]}, index=[i])
        res = pd.concat([res, result_df])
    adj = multipletests(res['pvalue'],alpha=0.05,method='fdr_bh')[1]
    res.insert(loc=2,column='Adj.P',value=adj)
    res = res[['pvalue','Adj.P']]
    return res
result = pv(gse119)
print(result)