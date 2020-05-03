import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys


sys.argv = ['diffX_cluster.py']
XL_namebase = r"C:\Users\rolep\Documents\Naithani Lab\SDRLK_Expression_Data\Abiotic Data"
metrix = 'euclidean'
if len(sys.argv) > 1:
    if sys.argv[1] == '-m':
        metrix = sys.argv[2]
        if len(sys.argv) > 3:
            XL_Name = sys.argv[3]
    else:
        XL_Name = sys.argv[1]

XL_handle = open((XL_namebase + r"\E-GEOD-38023-A-AFFY-126-query-results-C.txt"), "r")
df_Diff_XL = ["","","",""]
df_Diff_XL[0] = pd.read_csv(XL_handle,
                            index_col=0,
                            na_values="",
                            )
XL_handle.close()
df_Diff_XL[0].drop(labels=["Gene Name", "Design Element"], axis=1, inplace=True)
XL_handle = open((XL_namebase + r"\E-GEOD-41647-A-AFFY-126-query-results-C.txt"), "r")
df_Diff_XL[1] = pd.read_csv(XL_handle,
                            index_col=0,
                            na_values="",
                            )
XL_handle.close()
df_Diff_XL[1].drop(labels=["Gene Name", "Design Element"], axis=1, inplace=True)
XL_handle = open((XL_namebase + r"\E-MTAB-4994-A-AFFY-126-query-results-C.txt"), "r")
df_Diff_XL[2] = pd.read_csv(XL_handle,
                            index_col=0,
                            na_values="",
                            )
XL_handle.close()
df_Diff_XL[2].drop(labels=["Gene Name", "Design Element"], axis=1, inplace=True)
XL_handle = open((XL_namebase + r"\E-MTAB-5941-query-results-C.txt"), "r")
df_Diff_XL[3] = pd.read_csv(XL_handle,
                            index_col=0,
                            na_values="",
                            )
XL_handle.close()
df_Diff_XL[3].drop(labels=["Gene Name"], axis=1, inplace=True)

# df_Diff_XL.replace(to_replace="", value="0", inplace=True)
for t in range(0, 4):
    df_Diff_XL[t].fillna(value=0, inplace=True)
    df_cur = df_Diff_XL[t].filter(like='.foldChange')
    if t == 1:
        name = "E-GEOD-41647-new"
        curfig = [2.5, 5]
        metrica=metrix
        xticks=["B1", "B2", 'B3', 'B4']
    elif t == 2:
        name = "E-MTAB-4994-new"
        curfig = [1.5, 4]
        metrica=metrix
        xticks=['C1']
    elif t == 3:
        name = "E-MTAB-5941-new"
        curfig = [2.5, 8]
        metrica=metrix
        xticks=['D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7']
    else:
        name = "E-GEOD-38023-new"
        curfig = [2.5, 10]
        xticks=['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10']
        metrica=metrix

    Diff_XL_clus = sns.clustermap(df_cur,
                                  metric=metrica,
                                  cmap='RdBu',
                                  figsize=curfig,
                                  col_cluster=False,
                                  xticklabels=xticks,
                                  yticklabels=True,
                                  cbar_kws={"label": 'Log2 Expression Fold-Change',
                                            'orientation': 'horizontal'},
                                  center=0.0,
                                  # vmax=6.0
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
    Diff_XL_clus.ax_heatmap.tick_params(axis='y', labelsize=7)
    Diff_XL_clus.ax_heatmap.tick_params(axis='x', labelsize=10)
    plt.savefig(f"Diff_XL_SDRLK_abiotic_clus_{metrix}_{name}.png",
                bbox_inches='tight',
                )
