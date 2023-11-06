import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster, ward
data = pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\result_df.csv").rename(columns = {"0":"Index"})
cluster = pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\MCP_count_new_clu.csv")
col2keep=['Index','Ward_D']
cluster=cluster.loc[:,col2keep]
clean = data[['931','973']]
######
#compute classification score(PC1)
def cal_col4(row):
    return row['931'] * 0.74025 + row['973'] *  0.67233
clean['PC1'] = clean.apply(cal_col4, axis=1)
#compute classification score(PC2)
def cal_col5(row):
    return row['931'] * 0.672334895967201 + row['973'] *  -0.7402471125676703
clean['PC2'] = clean.apply(cal_col5, axis=1)
clean['Index']=data[['Index']]
merged_df = pd.merge(clean, cluster, on='Index')
col2keep=['PC1','PC2','Ward_D']
merged_df=merged_df.loc[:,col2keep]
# Perform PCA
fig = plt.figure(figsize=(6, 5.5))
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel('Principal Component 1(96.6%)', fontsize=15)
ax.set_ylabel('Principal Component 2(3.4%)', fontsize=15)
targets = [1,2]
colors = ['#F8766D','#00BFC4']
for target, color in zip(targets, colors):
    indicesToKeep = merged_df['Ward_D'] == target
    ax.scatter(merged_df.loc[indicesToKeep, 'PC1']
               , merged_df.loc[indicesToKeep, 'PC2']
               , c=color
               , s=40)
ax.legend(['Hot Tumor','Cold Tumor'])
plt.axvline(2.62333604530576,color='r')
ax.grid()
plt.show()

