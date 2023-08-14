import pandas as pd
File_dir=r"C:\Users\Heather P\Downloads\gdc_sample_sheet.2021-05-30.tsv"
sample_file = pd.read_csv(File_dir, sep='\t')
tumor_id = sample_file[sample_file['Sample Type'].str.contains('Primary Tumor')]
tumor_id.to_csv('tumor_sample',index=False)

################
# Sample DataFrame
data = {'1634': [1, 2, 3],
        '933': [4, 5, 6],
        '923': [7, 8, 9],
        '7075': [10, 11, 12]}
result_df = pd.DataFrame(data, index=['A', 'B', 'C'])

# List of column names
Tcells_id = ['933', '923']

# Calculate the sum and average for the specified columns in each row
for index in result_df.index:
    total_sum = 0
    num_columns = 0
    for col_name in result_df.columns:
        if col_name in Tcells_id:
            total_sum += result_df.loc[index, col_name]
            num_columns += 1
    if num_columns > 0:
        average_value = total_sum / num_columns
        print(f"Average value for columns in row {index}: {average_value}")
