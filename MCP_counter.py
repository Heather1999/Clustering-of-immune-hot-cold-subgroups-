import os
import pandas as pd
import numpy as np

csv_dir=r"C:\Users\Heather P\Desktop\github\T1\result_df.csv"
# Read the CSV file without header
result_df=pd.read_csv(csv_dir, header=None)
# Extract the first row as a list of strings
header_str = result_df.iloc[1, :].tolist()
# Set the first row as the column names
result_df.columns = header_str
# Remove the first row from the DataFrame
result_df= result_df.iloc[2:, :]
# Set the first column as the index
result_df.set_index(result_df.columns[0], inplace=True)

file_path = r"C:\Users\Heather P\Desktop\github\T1\109gene.csv"
gene_file=pd.read_csv(file_path)
###################################################################################################################
##pre counter classification(select all the genes belonging to specific cell type)
cell_pop=gene_file.iloc[:,2:4]
#1.Tcell_Entrezid
Tcells_id=[]
for i in range(len(cell_pop)):
    for j in range(len(cell_pop.columns)):
        if cell_pop.iloc[i,j] == "T cells":
            Tcells_id.append(cell_pop.iloc[i,j+1])
Tcells_id = list(map(str, Tcells_id))
#2.Cytotoxic lymphocytes
CL_id=[]
for i in range(len(cell_pop)):
    for j in range(len(cell_pop.columns)):
        if cell_pop.iloc[i,j] == "Cytotoxic lymphocytes":
            CL_id.append(cell_pop.iloc[i,j+1])
CL_id = list(map(str, CL_id))
#3.B lineage
BL_id=[]
for i in range(len(cell_pop)):
    for j in range(len(cell_pop.columns)):
        if cell_pop.iloc[i,j] == "B lineage":
            BL_id.append(cell_pop.iloc[i,j+1])
BL_id = list(map(str, BL_id))
#4.NK cells
NK_id=[]
for i in range(len(cell_pop)):
    for j in range(len(cell_pop.columns)):
        if cell_pop.iloc[i,j] == "NK cells":
            NK_id.append(cell_pop.iloc[i,j+1])
NK_id = list(map(str, NK_id))
#5.Monocytic lineage
ML_id=[]
for i in range(len(cell_pop)):
    for j in range(len(cell_pop.columns)):
        if cell_pop.iloc[i,j] == "Monocytic lineage":
            ML_id.append(cell_pop.iloc[i,j+1])
ML_id = list(map(str, ML_id))
#6.Myeloid dendritic cells
MD_id=[]
for i in range(len(cell_pop)):
    for j in range(len(cell_pop.columns)):
        if cell_pop.iloc[i,j] == "Myeloid dendritic cells":
            MD_id.append(cell_pop.iloc[i,j+1])
MD_id = list(map(str, MD_id))
#7.Neutrophils
NP_id=[]
for i in range(len(cell_pop)):
    for j in range(len(cell_pop.columns)):
        if cell_pop.iloc[i,j] == "Neutrophils":
            NP_id.append(cell_pop.iloc[i,j+1])
NP_id = list(map(str, NP_id))
#8.Endothelial cells
EC_id=[]
for i in range(len(cell_pop)):
    for j in range(len(cell_pop.columns)):
        if cell_pop.iloc[i,j] == "Endothelial cells":
            EC_id.append(cell_pop.iloc[i,j+1])
EC_id = list(map(str, EC_id))
#9.CD8 T cells
CD8_id=[]
for i in range(len(cell_pop)):
    for j in range(len(cell_pop.columns)):
        if cell_pop.iloc[i,j] == "CD8 T cells":
            CD8_id.append(cell_pop.iloc[i,j+1])
CD8_id = list(map(str, CD8_id))
#10.Fibroblasts
FB_id=[]
for i in range(len(cell_pop)):
    for j in range(len(cell_pop.columns)):
        if cell_pop.iloc[i,j] == "Fibroblasts":
            FB_id.append(cell_pop.iloc[i,j+1])
FB_id = list(map(str, FB_id))
###################################################################################################################
######MCP counter function
# Create a blank DataFrame to store the results
MCP_count_df = pd.DataFrame()
# Ensure column names are consistent
result_df.columns = result_df.columns.astype(str).str.split('.').str[0]
# Calculate the sum and average for the specified columns in each row
def average(a):
    global MCP_count_df  # Declare MCP_count_df as global
    for index in result_df.index:
        total_sum = 0
        num_columns = 0
        for col_name in result_df.columns:
            if col_name in a:
                total_sum += result_df.loc[index, col_name]
                num_columns += 1
        if num_columns > 0:
                average_value = total_sum / num_columns
                df2 = pd.DataFrame({'Index': [index], 'Cell Type': [col_name], 'Average': [average_value]})
                # Append the data to the MCP_count_df
                MCP_count_df = pd.concat([MCP_count_df, df2],ignore_index=True)
    return MCP_count_df
# Print the resulting DataFrame
average(Tcells_id)
average(CD8_id)
average(CL_id)
average(BL_id)
average(FB_id)
average(ML_id)
average(NK_id)
average(NP_id)
average(EC_id)
average(MD_id)
###Rearrange the data
cell_types = ['T cell', 'CD8', 'CL', 'BL', 'FB', 'ML', 'NK', 'NP', 'EC', 'MD']
start_indices = range(0, 360, 36)
for cell_type, start_index in zip(cell_types, start_indices):
    end_index = start_index + 36
    selected_data = MCP_count_df.iloc[start_index:end_index, 2].reset_index(drop=True)
    MCP_count_df[cell_type] = selected_data
# Drop columns using .drop() method
columns_to_drop = ['Cell Type', 'Average']  # List of column names to drop
rows_to_drop = range(36, 360)  # +1 to include the end_index
MCP_count_df = MCP_count_df.drop(columns=columns_to_drop)
MCP_count_df = MCP_count_df.drop(index=rows_to_drop)
print(MCP_count_df)