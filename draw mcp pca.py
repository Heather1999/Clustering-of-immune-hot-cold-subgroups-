import pandas as pd
import os
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
df = pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\MCP_count_new_clu.csv")
Type = list(df.columns[1:-1])

# Separating out the genes
x = df.loc[:,Type]
# Separating out the cluster
y = df.loc[:, ['Ward_D']].values
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
finalDf = pd.concat([finalDf, df['Ward_D'], df['Index']], axis=1)
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