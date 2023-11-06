import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import numpy as np
tcga_group = pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\MCP_count_new_clu.csv")
data = pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\result_df.csv").rename(columns = {"0":"Index"})
tcga_group = pd.merge(data,tcga_group, on='Index')
tcga_group = tcga_group[['Index','931','973','Ward_D']]
# Perform PCA
features = list(tcga_group.columns[1:3])
columns_to_keep=['931','973']
x=tcga_group.loc[:,columns_to_keep]
pca_2 = PCA(n_components=2)
# Separating out the features
x = tcga_group.loc[:, features]
# Separating out the target
y = tcga_group.loc[:, ['Ward_D']].values
x_fit = pca_2.fit(x)
eigen = pd.DataFrame(pca_2.components_)
path=r"C:\Users\Heather P\Desktop\github\T1\TCGA2geneeigen"
eigen.to_csv(path)
#======== TCGA 2 gens PCA
