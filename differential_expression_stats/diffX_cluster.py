import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

XL_Name = r"C:\Users\rolep\Documents\Naithani Lab\Biotic-SDRLK-data-for-fig-08-21-2019-Daemon-edited.csv"
sys.argv = ['diffX_cluster.py']
metrix = 'euclidean'

XL_handle = open(XL_Name, "r")
df_Diff_XL = pd.read_csv(XL_handle,
                         index_col=0,
                         na_values="",
                         )
XL_handle.close()

# df_Diff_XL.replace(to_replace="", value="0", inplace=True)
df_Diff_XL.fillna(value=0, inplace=True)

Diff_XL_clus = sns.clustermap(df_Diff_XL,
                              metric=metrix,
                              cmap='Spectral',
                              figsize=[8, 18],
                              col_cluster=False,
                              xticklabels=True,
                              yticklabels=True,
                              cbar_kws={"label": 'Log2 Expression Fold-Change',
                                        'orientation': 'horizontal'},
                              center=0.0
                              )
Diff_XL_clus.ax_col_dendrogram.set_visible(False)
DiffDendroBox = Diff_XL_clus.ax_col_dendrogram.get_position()
DiffDendroBox.y0 = (DiffDendroBox.y0 + 9 * DiffDendroBox.y1) / 10
DiffDendroWid = DiffDendroBox.y0 - DiffDendroBox.y1
DiffDendroBox.y0 = Diff_XL_clus.ax_col_dendrogram.get_position().y0
DiffDendroBox.y1 = DiffDendroBox.y0 - DiffDendroWid
Diff_XL_clus.cax.set_position(DiffDendroBox)
Diff_XL_clus.cax.xaxis.set_ticks_position("top")
Diff_XL_clus.cax.xaxis.set_label_position("top")
Diff_XL_clus.ax_heatmap.set_ylabel("Gene IDs")
plt.savefig(f"Diff_XL_SDRLK_clus_{metrix}.png",
            bbox_inches='tight',
            )
