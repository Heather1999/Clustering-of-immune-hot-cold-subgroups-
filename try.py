import pandas as pd
import numpy as np
file_dir=r"C:\Users\Heather P\Desktop\github\T1\result_df.csv"
min=pd.read_csv(file_dir)


A=[0.007252476,0.010896911,0.006900223,0.41193446,0.742771338,0.034275812]
def log2_function(value):
    return np.log2(value)
for i, each in enumerate(A):
    A[i] = log2_function(each)
A=sum(A)/6
B=0.00725247630646
B=log2_function(B)
print(B)