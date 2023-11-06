import pandas as pd
import os
import numpy as np
from sklearn.decomposition import PCA
pth=(r"C:\Users\Heather P\Desktop\github\T1\GSE119336_RNA.xlsx")
GSE119336=pd.read_excel(pth)
GSE119336.index=GSE119336.iloc[:,1]
GSE119336=GSE119336.iloc[:,3:]
pth=(r"C:\Users\Heather P\Desktop\github\T1\109gene.csv")
#gene=pd.read_csv(pth)
#gene=gene.loc[:,['Cell population','ENSEMBL ID']]
#row_to_filter = gene['ENSEMBL ID']
#GSE107943 = GSE107943.loc[GSE107943.index.isin(row_to_filter)]
for col in GSE119336.columns:
    min_value = float('inf')  # Initialize min_value with positive infinity
    for index, value in GSE119336[col].items():
        if index != 0:  # Exclude the first row
            numeric_value = float(value)
            if numeric_value != 0 and numeric_value < min_value:
                min_value = numeric_value 
    epsilon = 0.0000000001  # Adjust this threshold as needed
    for index, value in GSE119336[col].items():
        if index != 0:  # Exclude the first row
            numeric_value = float(value)
            if abs(numeric_value) < epsilon:
                GSE119336.at[index, col] = min_value
GSE119336.iloc[:, :] = GSE119336.iloc[:,:].applymap(lambda x: np.log2(x))

#GSE119336.to_csv(r"C:\Users\Heather P\Desktop\github\T1\GSE119336_TLS.csv")
