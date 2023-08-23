import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist
from sklearn import preprocessing
from scipy.cluster.hierarchy import linkage, fcluster, ward#############
from sklearn.decomposition import PCA
file_path = r"C:\Users\Heather P\Desktop\github\T1\MCP_count_new.csv"
MCP_score=pd.read_csv(file_path)
MCP_score.index=MCP_score.iloc[:,0]
MCP_score.columns = MCP_score.iloc[0,:]
MCP_score=MCP_score.iloc[1:,1:]
MCP_score = MCP_score.apply(pd.to_numeric)
file_path=r"C:\Users\Heather P\Desktop\github\T1\result_df.csv"
result_df=pd.read_csv(file_path)
result_df.index= result_df.iloc[:,0]
result_df=result_df.iloc[:,1:]
###########clustering then plotting
######
# Compute pairwise Euclidean distances using pdist
euclidean_dist = pdist(MCP_score.T, metric='euclidean')########
ward_d_linkage=linkage(euclidean_dist**(1/2),"ward")########
# Perform clustering based on a specified number of clusters
num_clusters = 2
MCP_score_cluster = MCP_score.T.copy()
MCP_score_cluster['Ward_D'] = fcluster(ward_d_linkage, num_clusters, criterion='maxclust')
file=r"C:\Users\Heather P\Desktop\github\T1\MCP_count_new_clu.csv"
MCP_score_cluster.to_csv(file)
###############
###############
#distance = pdist(MCP_score.values.T)
#linkage = ward(np.sqrt(distance))
#tumor_group = fcluster(linkage, t=2, criterion='maxclust')
#tumor_group = np.where(tumor_group==1, 1, 0)
#hot_group, cold_group = 1, 0
#col_colors = np.where(tumor_group==hot_group, "#FF8888", "#00BBCC")
#hot_indices, cold_indices = np.where(tumor_group == hot_group), np.where(tumor_group == cold_group)
#linkage[-1][0], linkage[-1][1] = linkage[-1][1], linkage[-1][0]
#linkage[-2][0], linkage[-2][1] = linkage[-2][1], linkage[-2][0]
#sns.clustermap(MCP_score, row_cluster=False, cmap="bwr", z_score=0, col_linkage=linkage, col_colors=col_colors)
#plt.show()
#############################PCA##109gene
scaled_data=preprocessing.scale(result_df)
pca=PCA()
pca.fit(scaled_data)#PCAmath place
pca_data=pca.transform(scaled_data)#generate coordinates ready to be plotted
per_var = np.round(pca.explained_variance_ratio_* 100, decimals=1)# calculate the percentage of each variation that PCA accounts for
labels = ['PC' + str(x) for x in range(1, len(per_var)+1)]
pca_df=pd.DataFrame(pca_data, columns=labels)
plt.scatter(pca_df.PC1,pca_df.PC2)
plt.xlabel('PC1 - {0}%'.format(per_var[0]))
plt.ylabel('PC2 - {0}%'.format(per_var[1]))
plt.show()
#############################PCA##MCP
#scaled_data=preprocessing.scale(MCP_score.T)
#pca=PCA()
#pca.fit(scaled_data)#PCAmath place
#pca_data=pca.transform(scaled_data)#generate coordinates ready to be plotted
##################create labels
#per_var = np.round(pca.explained_variance_ratio_* 100, decimals=1)# calculate the percentage of each variation that PCA accounts for
#labels = ['PC' + str(x) for x in range(1, len(per_var)+1)]
#############PCA plot
#pca_df=pd.DataFrame(pca_data, columns=labels)
#print(pca_df
#plt.scatter(pca_df.PC1,pca_df.PC2)
#plt.xlabel('PC1 - {0}%'.format(per_var[0]))
#plt.ylabel('PC2 - {0}%'.format(per_var[1]))
#plt.show()
