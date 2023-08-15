import os
import pandas as pd
import numpy as np
from collections import OrderedDict
folder_path = r"C:\Users\Heather P\Downloads\CHOL_raw\CHOL_raw"
dataframes = []  # To store DataFrames for each text file
csv_dir=r"C:\Users\Heather P\tumor_sample"
sample_filter=pd.read_csv(csv_dir)
####################open all files in one dataframe 
for filename in os.listdir(folder_path):
    if filename.endswith(".txt"):
        file_path = os.path.join(folder_path, filename)
        with open(file_path, "r") as file:
            content = file.read()
        # Extract the filename without extension
        filename_without_extension = filename.split('.')[0]
        df = pd.DataFrame({"File Name": [filename_without_extension], "Content": [content]})
        dataframes.append(df)
final_df = pd.concat(dataframes, ignore_index=True)
final_df.index=final_df.iloc[0:,0]
final_df=final_df.drop('File Name', axis=1)
###########################and change sample name
file_path=r"C:\Users\Heather P\Downloads\gdc_sample_sheet.2021-05-30.tsv"
sample_id=pd.read_csv(file_path,sep='\t')
for i in range(len(sample_id)):
    value = sample_id.iloc[i, 1]
    split_parts = value.split('.')
    sample_id.iloc[i, 1] = split_parts[0]
new_index = []
for value in final_df.index:
    # Find rows in sample_id where the 'File Name' column matches the value
    row_indices = sample_id.index[sample_id['File Name'] == value]   
    # Get the corresponding 'Sample ID' values from sample_id
    new_values = sample_id.iloc[row_indices, sample_id.columns.get_loc('Sample ID')].tolist()  
    # Append the values to the new_index list
    new_index.append(new_values)
# Create a new DataFrame from new_index
new_index_df = pd.DataFrame(new_index)
# Concatenate the new row and the original DataFrame
final_df.index=new_index_df.iloc[0:,0]
########################################################建立橫軸(基因)並且清乾淨.
# Split the content in the "Content" column into multiple columns
split_content = final_df["Content"].str.split("[\t\n]", expand=True)
# Get even columns and store gene names
even_columns = split_content.iloc[1, ::2]
gene_names = even_columns.drop(120966, axis=0).values.tolist()
gene_names = [item.split('.')[0] for item in gene_names]
# Remove even columns from split_content
split_content = split_content.drop(columns=even_columns.index)
split_content.columns=gene_names
final_df=split_content

####################################################篩選癌症sample(得到36)
#filtering tumor samples
file_id = sample_filter['Sample ID'].tolist()
final_df = final_df.loc[final_df.index.isin(file_id)]
################################################選出109基因
File_direc=r"C:\Users\Heather P\Downloads\genes.txt"
gene_file = pd.read_csv(File_direc, sep='\t')
valid_column_names = gene_file['ENSEMBL ID'].tolist()
valid_column_names = [item for item in valid_column_names if item != 'ENSG00000275037']
del valid_column_names[10]
final_df = final_df[valid_column_names]
####################置換基因呈現方式
# change ENSEMBL to ENTRENZID
gene = []
for value in final_df.columns:
    # Find rows in gene_file where the 'ENSEMBL ID' column matches the value
    row_indices = gene_file.index[gene_file['ENSEMBL ID'] == value]    
    # Get the corresponding 'ENTREZID' values from gene_file using .loc
    new_values = gene_file.loc[row_indices, 'ENTREZID'].tolist()
    # Append the values to the gene list
    gene.append(new_values)
gene = [values for values in gene if values]
gene = [item for sublist in gene for item in sublist]
gene = [str(item) for item in gene]
final_df.columns=gene
##################finding minimum value
# Iterate through each column to find the minimal value
for col in final_df.columns:
    min_value = float('inf')  # Initialize min_value with positive infinity
    for value in final_df[col]:
        numeric_value = float(value)
        if numeric_value != 0 and numeric_value < min_value:
            min_value = numeric_value  
    # Replace all zero values with the minimal value of the column
    epsilon = 0.0000000001  # Adjust this threshold as needed
    for i in range(len(final_df)):
        numeric_value = float(final_df.iloc[i][col])
        if abs(numeric_value) < epsilon:
            final_df.iloc[i, final_df.columns.get_loc(col)] = min_value
#####################log2
# Convert DataFrame to numeric type (ignore errors)
df_numeric = final_df.apply(pd.to_numeric)
# Define a function to calculate log2
def log2_function(value):
    return np.log2(value)
# Apply log2_function to each element of the DataFrame
final_df.iloc[:, :] = df_numeric.iloc[:, :].applymap(log2_function)
file_path = r"C:\Users\Heather P\Desktop\github\T1\result_df.csv"
final_df.to_csv(file_path)





