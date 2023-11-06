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
print(MCP_score)
###########clustering then plotting
######
# Compute pairwise Euclidean distances using pdist
#euclidean_dist = pdist(MCP_score.T, metric='euclidean')########
#ward_d_linkage=linkage(euclidean_dist**(1/2), method='ward')########
# Perform clustering based on a specified number of clusters
#MCP_score_a= fcluster(ward_d_linkage, 2, criterion='maxclust')
#hot_group, cold_group = 1, 2
#col_colors = np.where(tumor_group==hot_group, "#FF8888", "#00BBCC")###
#sns.clustermap(MCP_score, row_cluster=False, cmap="bwr", z_score=0,  col_linkage = MCP_score_a)
#plt.show()
##################################AUC ROC
##
from sklearn.metrics import roc_curve, auc
score = pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\PC1_score.csv")
tcga_group = pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\MCP_count_new_clu.csv")
tcga_group = tcga_group[['Index', 'Ward_D']]
score = score.merge(tcga_group, on="Index")
label = score['Ward_D'].replace({1: 1, 2: 0})
#把Ward_D變成機器學習可以理解的binary
###
#fpr, tpr, thresholds = roc_curve(label, score['PC1'])
#plt.figure()
#lw = 2
#roc_auc = auc(fpr, tpr)######
#plt.plot(fpr, tpr, color='darkorange',
#             lw=lw, label='ROC curve (area = %0.3f)' % roc_auc)
#plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
#plt.xlim([0.0, 1.0])
#plt.ylim([0.0, 1.05])
#plt.xlabel('False Positive Rate (FPR)')
#plt.ylabel('True Positive Rate (TPR)')
#plt.title('ROC Curve')
#plt.legend(loc='lower right')
#plt.show()
#####

def roc(PC):
    fpr, tpr, threshold = roc_curve(label, score[PC], drop_intermediate=False)
    roc_auc = auc(fpr, tpr)
    return fpr, tpr, roc_auc, threshold
# mark auc plot
auc_all = []
threshold_all = []
tpr_all = []
fpr_all = []

for i in list(score.columns[:-2]):
    fpr, tpr, roc_auc, threshold = roc(i)
    fpr_all.append(fpr)
    tpr_all.append(tpr)
    threshold_all.append(threshold)
    auc_all.append(roc_auc)
auc_all = pd.DataFrame(auc_all)
auc_all.rename(columns = {0:'AUC'},inplace=True)
no_gene = pd.Series(['gene'+str(i) for i in range(2,110)])
auc_all.insert(loc=0,column="NO. gene",value=no_gene)
# mark AUC plot (room in)
##
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
x = np.arange(2, len(auc_all) + 2)
ax.plot(x, auc_all['AUC'], label='AUC')
ax.legend(loc='lower right')
plt.ylabel('Area Under Curve')
plt.xlabel('Number of genes')
plt.xlim(0,20)
plt.ylim(0.9, )
plt.xticks(np.linspace(2, 20, num=10, dtype=int))
#plt.show()
# ##########AUC plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
x = np.arange(2, len(auc_all) + 2)
ax.plot(x, auc_all['AUC'], label='AUC')
ax.legend(loc='lower right')
plt.ylabel('Area Under Curve')
plt.xlabel('Number of genes')
plt.ylim(0,1.05)
plt.xticks(np.linspace(2, 109, num=10, dtype=int))
#plt.show()
##

