import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
dir_path = os.path.dirname(os.path.realpath(__file__))
ruta = os.path.join(dir_path,'heatmap_vals.csv')
df = pd.read_csv (ruta)
print(df)
df = df.round(decimals=2)
new_df = df.pivot(index="Value EB", columns="Value R", values="Genome diversity")
print(new_df)
ax = sns.heatmap(new_df)
plt.show()