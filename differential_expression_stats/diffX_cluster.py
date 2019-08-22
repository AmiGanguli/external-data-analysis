import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

XL_Name = r"C:\Users\rolep\Documents\Naithani Lab\SDRLK_Expression_Data\Biotic-data-file-02-for-Daemon-08-22-2019.csv"
sys.argv = ['diffX_cluster.py']
metrix = 'euclidean'

XL_handle = open(XL_Name, "r")
df_Diff_XL = ["","","","",""]
df_Diff_XL[0] = pd.read_csv(XL_handle,
                            index_col=0,
                            na_values="",
                            usecols=[0, 1, 2, 3],
                            skiprows=[0],
                            nrows=18,
)
XL_handle.close()
XL_handle = open(XL_Name, "r")
df_Diff_XL[1] = pd.read_csv(XL_handle,
                            index_col=0,
                            na_values="",
                            usecols=range(4, 23),
                            skiprows=[0],
                            nrows=43,
                            )
XL_handle.close()
XL_handle = open(XL_Name, "r")
df_Diff_XL[2] = pd.read_csv(XL_handle,
                            index_col=0,
                            na_values="",
                            usecols=range(23, 27),
                            skiprows=[0],
                            nrows=71,
                            )
XL_handle.close()
XL_handle = open(XL_Name, "r")
df_Diff_XL[3] = pd.read_csv(XL_handle,
                            index_col=0,
                            na_values="",
                            usecols=range(27, 31),
                            skiprows=[0],
                            nrows=41,
                            )
XL_handle.close()
XL_handle = open(XL_Name, "r")
df_Diff_XL[4] = pd.read_csv(XL_handle,
                            index_col=0,
                            na_values="",
                            usecols=range(31, 36),
                            skiprows=[0],
                            nrows=39,
                            )
XL_handle.close()

# df_Diff_XL.replace(to_replace="", value="0", inplace=True)
for t in range(0,5):
    df_Diff_XL[t].fillna(value=0, inplace=True)
    if t == 1:
        name = "X-oryzae-E-GEOD-36272"
        curfig = [8, 8]
    elif t == 2:
        name = "E-MTAB-4406"
        curfig = [3.5, 12]
    elif t == 3:
        name = "E-MTAB-5025"
        curfig = [3.5, 8]
    elif t == 4:
        name = "R. solani testing"
        curfig = [3.5, 8]
    else:
        name = "X1-E-GEOD-61833"
        curfig = [2.5, 4]

    Diff_XL_clus = sns.clustermap(df_Diff_XL[t],
                                  metric=metrix,
                                  cmap='RdBu',
                                  figsize=curfig,
                                  col_cluster=False,
                                  xticklabels=True,
                                  yticklabels=True,
                                  cbar_kws={"label": 'Log2 Expression Fold-Change',
                                            'orientation': 'horizontal'},
                                  center=0.0,
                                  #vmax=6.0
                                  )
    Diff_XL_clus.ax_col_dendrogram.set_visible(False)
    DiffDendroBox = Diff_XL_clus.ax_col_dendrogram.get_position()
    DiffDendroBox.y0 = (DiffDendroBox.y0 + 9 * DiffDendroBox.y1) / 10
    DiffDendroWid = DiffDendroBox.y0 - DiffDendroBox.y1
    DiffDendroBox.y0 = Diff_XL_clus.ax_col_dendrogram.get_position().y0
    DiffDendroBox.y1 = DiffDendroBox.y0 - 1.5*DiffDendroWid
    DiffDendroBox.y0 = DiffDendroBox.y1 + DiffDendroWid
    Diff_XL_clus.cax.set_position(DiffDendroBox)
    Diff_XL_clus.cax.xaxis.set_ticks_position("top")
    Diff_XL_clus.cax.xaxis.set_label_position("top")
    Diff_XL_clus.ax_heatmap.set_ylabel("Gene IDs")
    Diff_XL_clus.ax_heatmap.tick_params(labelsize=7)
    plt.savefig(f"Diff_XL_SDRLK_clus_{metrix}_{name}.png",
                bbox_inches='tight',
                )
