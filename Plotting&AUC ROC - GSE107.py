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
file_path = r"C:\Users\Heather P\Desktop\github\T1\MCP_count_GSE107943.csv"
MCP_score=pd.read_csv(file_path)
MCP_score.index=MCP_score['Index']
MCP_score=MCP_score.iloc[:,2:]
###########clustering then plotting
######
# Compute pairwise Euclidean distances using pdist
euclidean_dist = pdist(MCP_score, metric='euclidean')########
ward_d_linkage=linkage(euclidean_dist**(1/2), method='ward')########
# Perform clustering based on a specified number of clusters
tumor_group= fcluster(ward_d_linkage, 3, criterion='maxclust')
hot_group, cold_group, intermediate_group = 1, 2, 3
# Assign colors to the clusters
cluster_colors = {
    hot_group: "#00BBCC",
    cold_group: "#FF8888",
    intermediate_group: "#FDFEBA"}
col_colors = np.array([cluster_colors[label] for label in tumor_group])
#sns.clustermap(MCP_score.T, row_cluster=False, cmap="bwr", z_score=0, col_linkage=ward_d_linkage, col_colors=col_colors)
#plt.show()
MCP_score['group']=tumor_group
##############################PCA
#Type = list(MCP_score.columns[:-1])
# Separating out the genes
#x = MCP_score.loc[:,Type]
# Separating out the cluster
#y = MCP_score.loc[:, ['group']].values
#pca = PCA(n_components=2)
#x_fit = pca.fit(x)
#eigen = pd.DataFrame(pca.components_)
#PC1 = 0
#PC2 = 0
#for i in range(len(x.columns)):
#   PC1 = PC1 + (pd.Series(x.iloc[:, i] * eigen.iloc[0, i]))#每個sample乘上自己的eigenvector
#    PC2 = PC2 + (pd.Series(x.iloc[:, i] * eigen.iloc[1, i]))
#frame = {"PC1":PC1,"PC2":PC2}
#finalDf = pd.DataFrame(frame)
#finalDf = pd.concat([finalDf, MCP_score['group']], axis=1)
#finalDf['PC1'] = -finalDf['PC1']
#finalDf['PC2'] = -finalDf['PC2']
##2D
#fig = plt.figure(figsize=(6, 5.5))
#ax = fig.add_subplot(1, 1, 1)
#ax.set_xlabel('Principal Component 1 (%.1f'%(pca.explained_variance_ratio_[0]*100)+'%)', fontsize=15)
#ax.set_ylabel('Principal Component 2 (%.1f'%(pca.explained_variance_ratio_[1]*100)+'%)', fontsize=15)
#targets = [1,2,3]
#colors = ['#00BFC4','#F8766D','#FEC670']
#for target, color in zip(targets, colors):
#    indicesToKeep = finalDf['group'] == target
#    ax.scatter(finalDf.loc[indicesToKeep, 'PC1']
#               , finalDf.loc[indicesToKeep, 'PC2']
#               , c=color, s=40)
#ax.legend(['Hot Tumor','Cold Tumor','Intermediate'])
#ax.grid()
#ax.get_yaxis().set_visible(False)
#ax.invert_yaxis()
#plt.axvline(optimal_threshold['gene109'][0],color='r')#用roc算出來的optimal threshold
#plt.show()
###two gene model for this dataset
data = pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\GSE107943_109.csv")
data.index=data.iloc[:,0]
data=data.loc[['ENSG00000105369','ENSG00000156738'],:]
data=data.T
data=data.iloc[1:,:]
data.index.name='Index'
b = pd.merge(data,MCP_score, on='Index')
b = b[['ENSG00000105369','ENSG00000156738','group']]

# Perform PCA
features = ['ENSG00000105369','ENSG00000156738']
x = b.loc[:, features]
def cal_col4(row):
    return row['ENSG00000105369'] * 0.74025 + row['ENSG00000156738'] *  0.67233
x['PC1'] = x.apply(cal_col4, axis=1)
#compute classification score(PC2)
def cal_col5(row):
    return row['ENSG00000105369'] * 0.672334895967201 + row['ENSG00000156738'] *  -0.7402471125676703
x['PC2'] = x.apply(cal_col5, axis=1)
b = pd.merge(b,x, on='Index')
b=b.loc[:,['group','PC1','PC2']]
# Perform PCA
plt.figure(figsize=(10, 1))  # Set the figure size as needed, with a very small height
plt.xlabel('Principal Component 1', fontsize=15)
y_values = np.zeros(len(b))
legend_handles = [
    Line2D([0], [0], marker='o', color='w', label='Hot Tumor', markerfacecolor='red', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Cold Tumor', markerfacecolor='#00BFC4', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Intermediate', markerfacecolor='#FEC670', markersize=8)]
colors = ['#F8766D' if group == 2 else ('#00BFC4' if group == 1 else '#FEC670') for group in b['group']]
plt.scatter(b['PC1'], y_values, c=colors, s=40)
plt.axvline(1.231, color='r')
plt.legend(handles=legend_handles)
plt.grid()
#plt.show()
