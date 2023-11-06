import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import numpy as np
data = pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\min.csv",header = None)
data.columns=data.iloc[0,:]
data.iloc[:,1:] = np.log2(data.iloc[:,1:].astype('float'))
data = data.iloc[1:,:]
data.columns.values[0] = "Index"
tcga_group = pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\MCP_count_new_clu.csv")
tcga_group = tcga_group[['Index','Ward_D']]
tcga_group = pd.merge(data,tcga_group, on='Index')
df = tcga_group.drop(['Index'], axis=1)
# Perform PCA
x = df.drop(columns=['Ward_D'])  # Drop the target column
pca = PCA(n_components=2)
#x_pca = pca.fit_transform(x)
pca_109 = PCA(n_components=2).fit(x)
pc1_109, pc2_109 = np.dot(pca_109.components_, x.T)
# Create a DataFrame with PCA results
pca_df = pd.DataFrame()
pca_df['PC1']=pc1_109
pca_df['PC2']=pc2_109
pca_df = pd.concat([pca_df, tcga_group[['Ward_D', 'Index']]], axis=1)
#file = r"C:\Users\Heather P\Desktop\github\T1\pca_score.csv"
#pca_df.to_csv(file)
# Visualize PCA results
fig, ax = plt.subplots(figsize=(6, 5.5))
ax.set_xlabel(f'Principal Component 1 ({pca_109.explained_variance_ratio_[0]*100:.1f}%)', fontsize=15)
ax.set_ylabel(f'Principal Component 2 ({pca_109.explained_variance_ratio_[1]*100:.1f}%)', fontsize=15)

targets = [1, 2]
colors = ['#F8766D', '#00BFC4']
for target, color in zip(targets, colors):
    indicesToKeep = pca_df['Ward_D'] == target
    ax.scatter(pca_df.loc[indicesToKeep, 'PC1'], pca_df.loc[indicesToKeep, 'PC2'], c=color, s=40)

ax.legend(['Hot Tumor', 'Cold Tumor'])
ax.grid()
plt.axvline(1.06592,color='r')
plt.show()
# Calculate gene ranking
rotation = pd.DataFrame(pca_109.components_.T, columns=['PC1', 'PC2'])
rotation['Genesymbol'] = df.columns[:-1]
rotation['PC1'] = abs(rotation['PC1'])
gene_ranking = rotation.sort_values(by='PC1', ascending=False)

# Loop for gene ranking scores
score = pd.DataFrame()
for i in range(2, 110):
    data = tcga_group.loc[:, gene_ranking.iloc[0:i, -1]]
    pca = PCA(n_components=2)
    x_fit_transform = pd.DataFrame(pca.fit_transform(data), columns=['gene' + str(i), 'PC2'])
    eigen = pd.DataFrame(pca.components_)
    PC1 = (data * eigen.iloc[0]).sum(axis=1)
    score = pd.concat([score, pd.DataFrame(PC1, columns=['gene' + str(i)])], axis=1)
score = pd.concat([score, tcga_group[['Index']]], axis=1)
# Save and read scores
#score.to_csv(r"C:\Users\Heather P\Desktop\github\T1\PC1_score.csv", index=False)


