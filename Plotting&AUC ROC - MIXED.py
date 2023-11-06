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
euclidean_dist = pdist(MCP_score, metric='euclidean')########
ward_d_linkage=linkage(euclidean_dist**(1/2), method='ward')########
ward_d_linkage[-1][0], ward_d_linkage[-1][1] = ward_d_linkage[-1][1], ward_d_linkage[-1][0]
ward_d_linkage[-2][0], ward_d_linkage[-2][1] = ward_d_linkage[-2][1], ward_d_linkage[-2][0]
# Perform clustering based on a specified number of clusters
tumor_group= fcluster(ward_d_linkage, 2, criterion='maxclust')
tumor_group = np.where(tumor_group==1, 1, 0)
MCP_score['group']=tumor_group
##############################PCA
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
#legend_handles = [
#    Line2D([0], [0], marker='o', color='w', label='Hot Tumor', markerfacecolor='#E18683', markersize=8),
#    Line2D([0], [0], marker='o', color='w', label='Cold Tumor', markerfacecolor='#5581B0', markersize=8)]
#colors = ['#E18683' if group == 0 else '#5581B0' for group in b['group']]
####################################################107
###two gene model for this dataset
data = pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\GSE107943_109.csv")
MCP_score=pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\MCP_count_GSE107943.csv")
data.index=data.iloc[:,0]
data=data.loc[['ENSG00000105369','ENSG00000156738'],:]
data=data.T
data=data.iloc[1:,:]
data.index.name='Index'
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
MCP_score['group']=tumor_group
c= pd.merge(data,MCP_score, on='Index')
c = c[['ENSG00000105369','ENSG00000156738','group']]
# Perform PCA
features = ['ENSG00000105369','ENSG00000156738']
x = c.loc[:, features]
def cal_col4(row):
    return row['ENSG00000105369'] * 0.74025 + row['ENSG00000156738'] *  0.67233
x['PC1'] = x.apply(cal_col4, axis=1)
#compute classification score(PC2)
def cal_col5(row):
    return row['ENSG00000105369'] * 0.672334895967201 + row['ENSG00000156738'] *  -0.7402471125676703
x['PC2'] = x.apply(cal_col5, axis=1)
c = pd.merge(c,x, on='Index')
c=c.loc[:,['group','PC1','PC2']]
c.loc[c.iloc[:, 0] == 1, 0] = 2
c.loc[c.iloc[:, 0] == 2, 0] = 4
c.loc[c.iloc[:, 0] == 3, 0] = 3
c.columns=['group','PC1','PC2','group']
c=c.iloc[:,1:]
combined_data = pd.concat([b, c], ignore_index=True)
# Create the plot
plt.figure(figsize=(10, 1))  # Set the figure size as needed, with a very small height
plt.xlabel('Principal Component 1', fontsize=15)
y_values = np.zeros(len(combined_data))

# Define colors or markers for each part
color_mapping = {
    0: '#E18683',
    1: '#5581B0',
    4: '#A4312A',
    3: '#56BCF9',
    2: '#56BCF9',
    # Add more mappings as needed
}
colors = [color_mapping.get(group,'#56BCF9') for group in combined_data['group']]
marker_size = [40 if group == 0 else 40 for group in combined_data['group']]  # Adjust marker size as needed

# Plot the data for both parts on the same plot
plt.scatter(combined_data['PC1'], y_values, c=colors, s=marker_size)
plt.axvline(1.231, color='r')

# Add a legend
legend_handles = [
    Line2D([0], [0], marker='o', color='w', label='GSE119336-Hot', markerfacecolor='#E18683', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='GSE119336-Cold', markerfacecolor='#5581B0', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='GSE107943-Hot', markerfacecolor='#A4312A', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='GSE107943-Cold + Intermediate', markerfacecolor='#56BCF9', markersize=8)]
plt.legend(handles=legend_handles)

plt.grid()
plt.show()
