import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn
from scipy.cluster.hierarchy import dendrogram, linkage
file_path = r"C:\Users\Heather P\Desktop\github\T1\MCP_scores.csv"
MCP_score=pd.read_csv(file_path)
MCP_score.index=MCP_score.iloc[:,1]
columns_to_drop = MCP_score.iloc[:, 0:2]
MCP_score = MCP_score.drop(columns=columns_to_drop)
MCP_score = MCP_score.iloc[:, [0,1,2,6,3,5,9,7,8,4]]
MCP_score=MCP_score.T
print(MCP_score)
###########clustering then plotting
# Perform hierarchical clustering on columns
sn.clustermap(MCP_score,row_cluster=False, col_cluster=True,method='ward', z_score=0 ,metric='euclidean',cmap='bwr')
plt.show()