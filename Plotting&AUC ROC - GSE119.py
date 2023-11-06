import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
from scipy.spatial.distance import pdist
from sklearn import preprocessing
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, fcluster, ward#############
from sklearn.decomposition import PCA
file_path = r"C:\Users\Heather P\Desktop\github\T1\MCP_count_GSE119.csv"
MCP_score=pd.read_csv(file_path)
MCP_score.index=MCP_score['Index']
MCP_score=MCP_score.iloc[:,2:]
###########clustering then plotting
######
# Compute pairwise Euclidean distances using pdist
#euclidean_dist = pdist(MCP_score, metric='euclidean')########
#ward_d_linkage=linkage(euclidean_dist**(1/2), method='ward')########
#ward_d_linkage[-1][0], ward_d_linkage[-1][1] = ward_d_linkage[-1][1], ward_d_linkage[-1][0]
#ward_d_linkage[-2][0], ward_d_linkage[-2][1] = ward_d_linkage[-2][1], ward_d_linkage[-2][0]
# Perform clustering based on a specified number of clusters
#tumor_group= fcluster(ward_d_linkage, 2, criterion='maxclust')
MCP_score = MCP_score.sort_values(by='2-gene group')
tumor_group = np.where(MCP_score['2-gene group']==1, 1, 0)
hot_group, cold_group = 0, 1
col_colors = np.where(tumor_group==hot_group, "#FF8888", "#00BBCC")###
hot_indices, cold_indices = np.where(tumor_group == hot_group), np.where(tumor_group == cold_group)
sns.clustermap(MCP_score.T, row_cluster=False, cmap="bwr", z_score=0,  col_cluster=False,col_colors=col_colors)
plt.show()
MCP_score['group']=tumor_group
print(MCP_score)
##############################PCA
Type = list(MCP_score.columns[:-1])
# Separating out the genes
x = MCP_score.loc[:,Type]
# Separating out the cluster
y = MCP_score.loc[:, ['group']].values
pca = PCA(n_components=2)
x_fit = pca.fit(x)
eigen = pd.DataFrame(pca.components_)
PC1 = 0
PC2 = 0
for i in range(len(x.columns)):
    PC1 = PC1 + (pd.Series(x.iloc[:, i] * eigen.iloc[0, i]))#每個sample乘上自己的eigenvector
    PC2 = PC2 + (pd.Series(x.iloc[:, i] * eigen.iloc[1, i]))
frame = {"PC1":PC1,"PC2":PC2}
finalDf = pd.DataFrame(frame)
finalDf = pd.concat([finalDf, MCP_score['group']], axis=1)
finalDf['PC1'] = -finalDf['PC1']
finalDf['PC2'] = -finalDf['PC2']
##2D
#fig = plt.figure(figsize=(6, 5.5))
#ax = fig.add_subplot(1, 1, 1)
#ax.set_xlabel('Principal Component 1 (%.1f'%(pca.explained_variance_ratio_[0]*100)+'%)', fontsize=15)
##ax.set_ylabel('Principal Component 2 (%.1f'%(pca.explained_variance_ratio_[1]*100)+'%)', fontsize=15)
##targets = [0,1]
#colors = ['#F8766D','#00BFC4']
##for target, color in zip(targets, colors):
#    indicesToKeep = finalDf['group'] == target
#    ax.scatter(finalDf.loc[indicesToKeep, 'PC1']
#               , finalDf.loc[indicesToKeep, 'PC2']
#               , c=color, s=40)
#ax.legend(['Hot Tumor','Cold Tumor'])
#ax.grid()
#ax.get_yaxis().set_visible(False)
#ax.invert_yaxis()
#plt.axvline(optimal_threshold['gene109'][0],color='r')#用roc算出來的optimal threshold
###two gene model for this dataset
data = pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\GSE119336_105RNA.csv")
data.index=data.iloc[:,0]
data=data.loc[['MS4A1','CD79A'],:]
data=data.T
data = data.rename_axis("Index")
b = pd.merge(data,MCP_score, on='Index')
b = b[['MS4A1','CD79A','group']]
# Perform PCA
features = ['MS4A1','CD79A']
x = b.loc[:, features]
def cal_col4(row):
    return row['MS4A1'] * 0.74025 + row['CD79A'] *  0.67233
x['PC1'] = x.apply(cal_col4, axis=1)
#compute classification score(PC2)
def cal_col5(row):
    return row['MS4A1'] * 0.672334895967201 + row['CD79A'] *  -0.7402471125676703
x['PC2'] = x.apply(cal_col5, axis=1)
b = pd.merge(b,x, on='Index')
b=b.loc[:,['group','PC1','PC2']]
print(b)
# Perform PCA
plt.figure(figsize=(10, 1))  # Set the figure size as needed, with a very small height
plt.xlabel('Principal Component 1', fontsize=15)
y_values = np.zeros(len(b))
legend_handles = [
    Line2D([0], [0], marker='o', color='w', label='Hot Tumor', markerfacecolor='red', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Cold Tumor', markerfacecolor='#00BFC4', markersize=8)]
colors = ['#F8766D' if group == 0 else '#00BFC4' for group in b['group']]
plt.scatter(b['PC1'], y_values, c=colors, s=40)
plt.axvline(1.231, color='r')
plt.legend(handles=legend_handles)
plt.grid()
plt.show()
