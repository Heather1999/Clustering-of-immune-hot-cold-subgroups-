import pandas as pd
import os
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
#optimal_threshold = pd.read_csv(r"C:\Users\Heather P\Downloads\New folder1\cheat\optimal_threshold.csv")
data = pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\min.csv",header = None)
data.columns=data.iloc[0,:]
data.iloc[:,1:] = np.log2(data.iloc[:,1:].astype('float'))
data = data.iloc[1:,:]
data.columns.values[0] = "Index"

tcga_group = pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\MCP_count_new_clu.csv")
tcga_group = tcga_group[['Index','Ward_D']]
tcga_group = pd.merge(data,tcga_group, on='Index')

gene = pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\genes.txt",sep='\t',index_col=0)
df = tcga_group.drop(['Index'], axis=1)
genes = list(df.columns[:-1])

# Separating out the genes
x = df.loc[:, genes]
# Separating out the cluster
y = df.loc[:, ['Ward_D']].values
# Applying PCA to the gene data
pca = PCA(n_components=2)
x_fit = pca.fit(x)
eigen = pd.DataFrame(pca.components_)
PC1 = 0
PC2 = 0
for i in range(len(x.columns)):
    PC1 = PC1 + (pd.Series(x.iloc[:, i] * eigen.iloc[0, i]))#每個sample(gene?)乘上自己的eigenvector
    PC2 = PC2 + (pd.Series(x.iloc[:, i] * eigen.iloc[1, i]))
#pca.explained_variance_) 變異數 pca.mean_ # 成分平均 pca.explained_variance_ratio_ # 成分比例 pca.components_# 成分的投影軸向量，又稱特徵向量（eigenvectors）
frame = {"PC1":PC1,"PC2":PC2}
finalDf = pd.DataFrame(frame)

finalDf = pd.concat([finalDf, tcga_group['Ward_D'], tcga_group['Index']], axis=1)
##
fig = plt.figure(figsize=(6, 5.5))
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel('Principal Component 1 (%.1f'%(pca.explained_variance_ratio_[0]*100)+'%)', fontsize=15)
ax.set_ylabel('Principal Component 2 (%.1f'%(pca.explained_variance_ratio_[1]*100)+'%)', fontsize=15)
targets = [1,2]
colors = ['#F8766D','#00BFC4']
for target, color in zip(targets, colors):
    indicesToKeep = finalDf['Ward_D'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'PC1']
               , finalDf.loc[indicesToKeep, 'PC2']
               , c=color
               , s=40)
ax.legend(['Hot Tumor','Cold Tumor'])
ax.grid()
#plt.axvline(optimal_threshold['gene109'][0],color='r')#用roc算出來的optimal threshold
plt.show()


rotation = pca.components_  # = R pca的rotation
# pca.explained_variance_ratio_# = 每個pc解釋的比例
rotation = rotation.T

rotation = pd.concat([pd.DataFrame(rotation, columns=["PC1", "PC2"]), pd.DataFrame(df.columns[:-1], columns=['Genesymbol'])], axis=1)
rotation = rotation.iloc[:, [0, 2]]
rotation.iloc[:, 0] = abs(rotation.iloc[:, 0])
gene_ranking = rotation.sort_values(by='PC1', ascending=False)
#gene_ranking.reset_index(drop=True,inplace=True)
#gene_ranking.to_csv(r"C:\Users\Heather P\Desktop\github\T1\gene_ranking.csv",index=False)
##
def rank(gene_num):
    data = tcga_group.loc[:, list(gene_ranking.iloc[0:gene_num, 1])]
    data.reset_index(drop = True)
    pca = PCA(n_components=2)
    x_fit_transform = pd.DataFrame(pca.fit_transform(data), columns=['gene' + str(gene_num), 'PC2'])
    x_fit = pca.fit(data)
    eigen = pd.DataFrame(pca.components_)
    PC1 = 0
    PC2 = 0
    for i in range(len(data.columns)):
        PC1 = PC1 + (pd.Series(data.iloc[:, i] * eigen.iloc[0, i]))  # 每個sample乘上自己的eigenvector
        PC2 = PC2 + (pd.Series(data.iloc[:, i] * eigen.iloc[1, i]))
    return PC1

##
score = pd.DataFrame()
for i in range(2,110):
    score = pd.concat([score,pd.DataFrame(rank(i),columns=['gene'+str(i)])],axis=1)

score = pd.concat([score,tcga_group['Index']],axis=1)
score.to_csv(r"C:\Users\Heather P\Desktop\github\T1\PC1_score.csv",index=False)
