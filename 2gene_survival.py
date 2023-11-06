import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
gse=pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\GSE107943_survival - 2 gene.csv")
tcga=pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\TCGA survival_2gene - 複製.csv")
print(gse)
# Create separate DataFrames for different groups or subpopulations#1=high#931=MS4A1
gse['characteristics_ch1.5'] = gse['characteristics_ch1.5'].str.extract(r'(\d+)').astype(int)
gse['characteristics_ch1.3'] = gse['characteristics_ch1.3'].str.extract(r'(\d+)').astype(int)
group1_df = gse[gse['groupCD79A'] == 1]
group2_df = gse[gse['groupCD79A'] == 2]
#gse.loc[25,'group']=1
# Fit Kaplan-Meier curves for each group
group1_df = group1_df.dropna(subset=['dsfree(mo):ch1', 'characteristics_ch1.3'])
group2_df = group2_df.dropna(subset=['dsfree(mo):ch1', 'characteristics_ch1.3'])
kmf_group1 = KaplanMeierFitter()
kmf_group2 = KaplanMeierFitter()
kmf_group1.fit(durations=group1_df['dsfree(mo):ch1'], event_observed=group1_df['characteristics_ch1.3'], label='CD79A High')
kmf_group2.fit(durations=group2_df['dsfree(mo):ch1'], event_observed=group2_df['characteristics_ch1.3'], label='CD79A Low')
# Plot survival curves for both groups
ax = kmf_group1.plot(ci_show=False, color='#E7796F')
kmf_group2.plot(ax=ax,ci_show=False,color='#69C9CC')
plt.title('')
plt.xlabel('months')
plt.ylabel('Disease free interval prportion')
results = logrank_test(
    gse['dsfree(mo):ch1'][gse['groupCD79A'] == 1],
    gse['dsfree(mo):ch1'][gse['groupCD79A'] == 2],
    event_observed_A=gse['characteristics_ch1.3'][gse['groupCD79A'] == 1],
    event_observed_B=gse['characteristics_ch1.3'][gse['groupCD79A'] == 2]
)
print(results.p_value)
plt.show()