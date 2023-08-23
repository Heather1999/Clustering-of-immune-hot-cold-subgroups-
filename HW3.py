import pandas as pd
file = r"C:\Users\Heather P\Downloads\New folder1\.genes.results"
RSEM = pd.read_csv(file,delimiter="\t")
RSEM = RSEM [['gene_id','effective_length']]
file = r"C:\Users\Heather P\Downloads\New folder1\SRR3184285ReadsPerGene.out.tab"
STAR = pd.read_csv(file,delimiter='\t')
STAR = STAR.iloc[3:,0:2]
##############################Concatenate the two
RSEM['unstranded_count'] = ''
# Create a dictionary mapping 'N_unmapped' values to '376989' values
star_mapping = dict(zip(STAR['N_unmapped'], STAR['376989']))
# Use map to populate 'unstranded count' column
RSEM['unstranded_count'] = RSEM['gene_id'].map(star_mapping)
RSEM = RSEM.query('(unstranded_count.notna()) & (effective_length.notna()) & (unstranded_count != 0) & (effective_length != 0)')
#########################################################compute FPKM unstranded
sum_count = RSEM['unstranded_count'].sum()
#FPKM = unstranded count/(sum unstranded count)*0.000001*effective_length*0.001
RSEM['FPKM'] = RSEM['unstranded_count'] / (sum_count*0.000001*RSEM['effective_length']*0.001)
FPKM_sum = RSEM['FPKM'].sum()
#TPM
RSEM['TPM'] = RSEM['FPKM'] / FPKM_sum*1000000
file=r"C:\Users\Heather P\Downloads\New folder1\test.csv"
RSEM.to_csv(file,index=False)
print(RSEM.head())

