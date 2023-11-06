import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
from statannot import add_stat_annotation
gse107=pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\GSE107943_Box.csv")
gse119=pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\GSE119336_Box.csv")
tcga=pd.read_csv(r"C:\Users\Heather P\Desktop\github\T1\result_df - Box.csv")
#################################################################################tcga
#x = "Ward_D"
#y = "973"
#order = [1,2]
#ax = sns.boxplot(data=tcga, x=x, y=y, order=order,showfliers = False)
#add_stat_annotation(ax, data=tcga, x=x, y=y, order=order,
#                    box_pairs=[(1,2)],
#                    test='Mann-Whitney', text_format='simple', verbose=2)
#plt.ylim(0,20)
#plt.xlabel('Group')
#plt.ylabel(' ')
#plt.show()
##############################################################################gse
x = "group"
y = "CD79A"
gse119.index=gse119.iloc[:,0]
gse119=(gse119.iloc[:,1:]).T
print(gse119)
order = [1,2]
ax = sns.boxplot(data=gse119, x=x, y=y, order=order,showfliers = False)
add_stat_annotation(ax, data=gse119, x=x, y=y, order=order,
                    box_pairs=[(1,2)],
                    test='Mann-Whitney', text_format='simple', verbose=2)
plt.ylim(0,20)
plt.xlabel('group')
plt.ylabel(' ')
plt.show()