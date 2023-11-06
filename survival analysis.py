import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
gse = pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\GSE107943_survival.csv")
# Create separate DataFrames for different groups or subpopulations
gse['characteristics_ch1.5'] = gse['characteristics_ch1.5'].str.extract(r'(\d+)').astype(int)
gse['characteristics_ch1.3'] = gse['characteristics_ch1.3'].str.extract(r'(\d+)').astype(int)
gse.loc[gse['group'] == 3, 'group'] = 1
group1_df = gse[gse['group'] == 1]
group2_df = gse[gse['group'] == 2]
#gse.loc[25,'group']=1
# Fit Kaplan-Meier curves for each group
kmf_group1 = KaplanMeierFitter()
kmf_group2 = KaplanMeierFitter()
kmf_group1.fit(durations=group1_df['survival(mo):ch1'], event_observed=group1_df['characteristics_ch1.5'], label='Cold Tumor')
kmf_group2.fit(durations=group2_df['survival(mo):ch1'], event_observed=group2_df['characteristics_ch1.5'], label='Hot Tumor')
# Plot survival curves for both groups
ax = kmf_group1.plot(ci_show=False, color='#69C9CC')
kmf_group2.plot(ax=ax,ci_show=False,color='#E7796F')
plt.title('')
plt.xlabel('Days')
plt.ylabel('Disease free interval proportion')
results = logrank_test(
    gse['survival(mo):ch1'][gse['group'] == 1],
    gse['survival(mo):ch1'][gse['group'] == 2],
    event_observed_A=gse['characteristics_ch1.5'][gse['group'] == 1],
    event_observed_B=gse['characteristics_ch1.5'][gse['group'] == 2]
)
print(results.p_value)
plt.show()