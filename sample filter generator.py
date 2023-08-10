import pandas as pd
File_dir=r"C:\Users\Heather P\Downloads\gdc_sample_sheet.2021-05-30.tsv"
sample_file = pd.read_csv(File_dir, sep='\t')
tumor_id = sample_file[sample_file['Sample Type'].str.contains('Primary Tumor')]
print(sample_file)
tumor_id.to_csv('tumor_sample',index=False)

####站存
