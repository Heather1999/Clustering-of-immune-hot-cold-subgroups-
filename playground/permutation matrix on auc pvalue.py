import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
df={'hrsofstudying':[2,2,20,10,4,15,26],'grade':[0,0,1,1,0,1,1]}
df=pd.DataFrame(df)

observed_auc = roc_auc_score(df['grade'], df['hrsofstudying'])
# Number of permutations for the permutation test
num_permutations = 100
permuted_aucs = []
for _ in range(num_permutations):
    shuffled_label=np.random.permutation(df['grade'])
    shuffled_auc=roc_auc_score(shuffled_label,df['hrsofstudying'])
    permuted_aucs.append(shuffled_auc)
p_value=(np.sum(np.array(permuted_aucs) >= observed_auc)+1)/(num_permutations+1)

print(p_value)
