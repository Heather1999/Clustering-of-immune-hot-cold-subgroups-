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

# Get the first row as column names
#new_column_names = final_df.iloc[0]
# Rename columns
#final_df.columns = new_column_names
print(final_df)


# Transpose the DataFrame
#final_df = final_df.T

# Truncate everything after the decimal point in the index
#final_df.index = final_df.index.str.split('.').str[0]

#matching_mask = final_df.iloc[0].isin(sample_filter['File Name'].str.strip()).any()
