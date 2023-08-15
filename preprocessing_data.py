import os
import pandas as pd
import numpy as np
folder_path = r"C:\Users\Heather P\Downloads\CHOL_raw\CHOL_raw"
dataframes = []  # To store DataFrames for each text file
csv_dir=r"C:\Users\Heather P\tumor_sample"
sample_filter=pd.read_csv(csv_dir)

for filename in os.listdir(folder_path):
    if filename.endswith(".txt"):
        file_path = os.path.join(folder_path, filename)
        with open(file_path, "r") as file:
            content = file.read()
        df = pd.DataFrame({"File Name": [filename], "Content": [content]})
        dataframes.append(df)
final_df = pd.concat(dataframes, ignore_index=True)

# Split the content in the "Content" column into multiple columns
split_content = final_df["Content"].str.split("[\t\n]", expand=True)  # Split by both tabs and newlines
even_columns = split_content.iloc[:, ::2]
Gene_names = []
# Loop through each cell in the DataFrame and append values to the list
for column_name in even_columns.columns:
    Gene_names.extend(even_columns[column_name].tolist())
# Remove duplicates while preserving order
Gene_names = list(dict.fromkeys(Gene_names))
Gene_names = Gene_names[:-1]
# Drop even-numbered columns
even_columns = split_content.columns[::2]  # Select even-numbered columns using slicing
split_content = split_content.drop(columns=even_columns, axis=1)
# Add the list as the first row using .loc
split_content.loc[-1] = Gene_names
split_content.index = split_content.index + 1
split_content = split_content.sort_index()
# Combine the split content with the "File Name" column
final_df = pd.concat([final_df["File Name"], split_content], axis=1)
# Set the "File Name" column as the index
final_df.set_index("File Name", inplace=True)
index = final_df.index.tolist()
index.insert(0, index.pop())
final_df.index = index


#filtering tumor samples
file_id = sample_filter['File Name'].tolist()
modified_id = []
for file_id in file_id:
    parts = file_id.split(".")  # Split the string by "."
    desired_parts = parts[:-1]  # Remove the last part
    modified_string = ".".join(desired_parts)  # Join the parts back together using "."
    modified_id.append(modified_string)
# Filter the rows based on the filter_list
filtered_df = final_df.loc[final_df.index.isin(modified_id)]
#select my nan
nan_row=final_df.loc[np.nan]
nan_row=pd.DataFrame(nan_row)
# Add the NaN row to the top of the filtered DataFrame
result_df = pd.concat([nan_row.T, filtered_df], ignore_index=False)
#Get the first row and copy it
row_nan = result_df.iloc[0].copy()
# Loop through the cells in the first row and modify values
for column in row_nan.index:
    if isinstance(row_nan[column], str):
        row_nan[column] = row_nan[column].split('.')[0]
# Update the modified row back to the DataFrame
result_df.iloc[0] = row_nan


###########
File_direc=r"C:\Users\Heather P\Downloads\genes.txt"
gene_file = pd.read_csv(File_direc, sep='\t')
valid_column_names = gene_file['ENSEMBL ID'].tolist()
result_df=result_df.T
result_df = result_df.loc[result_df.iloc[:, 0].isin(valid_column_names)]
#109 gene name 
valid_column_names= [col_name for col_name in valid_column_names if col_name in result_df.iloc[:, 0].values]
gene_file = gene_file[gene_file['ENSEMBL ID'].isin(valid_column_names)]
############
# change ENSEMBL to ENTRENZID
result_df.iloc[:,0]=result_df.iloc[:,0].sort_values(ascending=True)
gene_file = gene_file.sort_values(by='ENSEMBL ID',ascending=True)
gene_xchange = gene_file.iloc[:, 2]
result_df.iloc[:,0]=gene_xchange
result_df=result_df.T
########################
result_df.index = result_df.index.str.split('.').str[0]
file_path=r"C:\Users\Heather P\Downloads\gdc_sample_sheet.2021-05-30.tsv"
sample_id=pd.read_csv(file_path,sep='\t')
for i in range(len(sample_id)):
    value = sample_id.iloc[i, 1]
    split_parts = value.split('.')
    sample_id.iloc[i, 1] = split_parts[0]
new_index = []
for value in result_df.index:
    # Find rows in sample_id where the 'File Name' column matches the value
    row_indices = sample_id.index[sample_id['File Name'] == value]   
    # Get the corresponding 'Sample ID' values from sample_id
    new_values = sample_id.iloc[row_indices, sample_id.columns.get_loc('Sample ID')].tolist()  
    # Append the values to the new_index list
    new_index.append(new_values)
# Create a new DataFrame from new_index
new_index_df = pd.DataFrame(new_index)
# Concatenate the new row and the original DataFrame
result_df.index=new_index_df.iloc[:,0]


#####finding minimum value
# Initialize the minimum value to a large number
min_value = float('inf')

# Iterate through each column to find the minimal value
for col in result_df.columns:
    min_value = float('inf')  # Initialize min_value with positive infinity
    for value in result_df[col]:
        numeric_value = float(value)
        if numeric_value != 0 and numeric_value < min_value:
            min_value = numeric_value  
    # Replace all zero values with the minimal value of the column
    epsilon = 0.0000000001  # Adjust this threshold as needed
    for i in range(len(result_df)):
        numeric_value = float(result_df.iloc[i][col])
        if abs(numeric_value) < epsilon:
            result_df.iloc[i, result_df.columns.get_loc(col)] = min_value


#####log2
# Convert DataFrame to numeric type (ignore errors)
df_numeric = result_df.apply(pd.to_numeric)
# Define a function to calculate log2
def log2_function(value):
    return np.log2(value)
# Apply log2_function to each element of the DataFrame
result_df.iloc[1:, :] = df_numeric.iloc[1:, :].applymap(log2_function)

file_path = r"C:\Users\Heather P\Desktop\github\T1\result_df.csv"
result_df.to_csv(file_path)
#file_path = r"C:\Users\Heather P\Desktop\github\T1\109gene.csv"
#gene_file.to_csv(file_path)

