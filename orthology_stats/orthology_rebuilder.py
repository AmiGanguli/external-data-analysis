#!C:\Program Files (x86)\Python
import json
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns
import copy

url_base = 'https://plantreactomedev.gramene.org/ContentService'
headers = {'accept': 'application/json'}

# customizable species list - r63
species_list = [# Monocots
                "Oryza sativa Indica", "Oryza sativa aus kasalath", "Oryza australiensis", "Oryza barthii", "Oryza brachyantha", "Oryza glaberrima", "Oryza glumaepatula", 
                "Oryza granulata", "Oryza longistaminata", "Oryza meridionalis", "Oryza minuta", "Oryza nivara", "Oryza officinalis", "Oryza punctata", "Oryza rufipogon", "Hordeum vulgare",
                "Setaria italica", "Sorghum bicolor", "Triticum aestivum", "Triticum turgidum", "Triticum urartu", "Triticum dicoccoides", "Zea mays", "Brachypodium distachyon",
                "Aegilops tauschii", "Eragrostis curvula", "Eragrostis tef", "Phyllostachys edulis", "Saccharum spontaneum", "Leersia perrieri",
                "Panicum hallii FIL2", "Panicum hallii var. hallii HAL2", "Zoysia japonica", "Ananas comosus", "Asparagus officinalis", "Musa acuminata",
                "Dioscorea rotundata", 
                # Dicots
                "Beta vulgaris",  "Eucalyptus grandis", "Coffea canephora", "Phoenix dactylifera", "Vitis vinifera", "Nelumbo nucifera", 
                "Helianthus annuus", "Cynara cardunculus", "Salvia hispanica", "Erythranthe guttata", "Olea europaea var. sylvestris", 
                "Ipomoea triloba", "Nicotiana attenuata", "Solanum tuberosum", "Solanum lycopersicum",  "Capsicum annuum", "Daucus carota", 
                "Arachis ipaensis", "Arachis duranensis", "Trifolium pratense", "Vigna radiata", "Vigna angularis", "Phaseolus vulgaris", 
                "Medicago truncatula", "Lupinus angustifolius", "Glycine max", "Cicer arietinum", "Cajanus cajan", "Capsella rubella", 
                "Brassica oleracea", "Brassica rapa", "Brassica napus", "Arabidopsis thaliana", "Arabidopsis lyrata", "Arabidopsis halleri", 
                "Populus trichocarpa", "Manihot esculenta", "Jatropha curcas", "Cucumis sativus", "Citrullus lanatus", "Corchorus capsularis",
                "Corchorus olitorius", "Gossypium raimondii", "Theobroma cacao", "Actinidia chinensis", "Fragaria vesca", "Prunus persica", 
                "Prunus avium", "Malus domestica", "Humulus lupulus primary", "Humulus lupulus haplotig", "Cannabis sativa", 
                "Cannabis sativa subsp. indica", "Pistacia vera", "Citrus clementina", "Citrus sinensis",
                # Basal angiosperm
                "Amborella trichopoda",
                # Gymnosperms
                "Pinus taeda", "Picea abies",
                # Lower Plants
                "Marchantia polymorpha", "Physcomitrella patens", "Chara braunii", "Selaginella moellendorffii", "Chlamydomonas reinhardtii", 
                 "Ostreococcus lucimarinus", "Galdieria sulphuraria", "Cyanidioschyzon merolae", "Chondrus crispus",
                "Synechocystis sp. PCC 6803",
                ]

# example customizations
#species_list = ["Arabidopsis thaliana"] # test run 0 - AT only - JP 20201224
#species_list = ["Arabidopsis thaliana","Musa acuminata","Triticum aestivum","Pinus taeda","Selaginella moellendorffii"] # test run 1 - multiclade - JP 20201224
#species_list = ["Selaginella moellendorffii","Oryza australiensis", "Oryza meyeriana var. granulata", "Oryza minuta", "Oryza officinalis"] # test run 2 - outlier - JP 20201224

def recurs_get_paths(sub_dict, path_list, term_path_list):
    pathTypes = ['Pathway', 'TopLevelPathway']
    rxnTypes = ['Reaction', 'BlackBoxEvent']
    if sub_dict['stId'] in path_list and sub_dict['type'] in pathTypes:
        termFlag = False
        for child in sub_dict['children']:
            if child['type'] == 'Pathway':
                path_list.append(child['stId'])
                # IDs[f"{child['stId']}"] = child['name']

                recurs_get_paths(child, path_list, term_path_list)
            elif child['type'] in rxnTypes and termFlag is False:
                termFlag = True
                sub_term_list = ["", ""]
                sub_term_list[0] = sub_dict['stId']
                sub_term_list[1] = sub_dict['name']
                term_path_list.append(sub_term_list)
                # IDs[f"{child['name']}"] = child['stId']
    elif sub_dict['type'] in pathTypes:
        for child in sub_dict['children']:
            if child['type'] in pathTypes:
                recurs_get_paths(child, path_list, term_path_list)
    return


# if __name__ == '__main__':
CountFlag = True

# something something interface for user to decide which pathways they want?
sys.argv = ['orthology_rebuilder.py', 'W', '-r', '-log',
            'R-OSA-2744345']
# P -> Protein Level; R -> Reaction Level; W -> pathWay Level; A -> All three
if sys.argv[1] not in ['P', 'R', 'W', 'A']:
    print("Sorry, that is not an option for group sorting level. The options are:\n"
          "\'P\': sort by Protein ID.\n"
          "\'R\': sort by Reaction ID.\n"
          "\'W\': sort by pathWay ID.\n"
          "\'A\': sort by All of the above.")
    exit()
levelFlag = sys.argv[1]
ratioFlag = False
if sys.argv[2] is '-r':
    ratioFlag = True
logFlag = False
if sys.argv[3] is '-log':
    logFlag = True
path_list = sys.argv[4:]
term_path_list = []

path_url = '/data/eventsHierarchy/4530'
full_url = url_base + path_url
response = requests.get(full_url, headers=headers).json()

base_dict = {}
for TopLevel in response:
    if TopLevel['name'] == "Metabolism and regulation":
        base_dict = TopLevel
        break
recurs_get_paths(base_dict, path_list, term_path_list)

# take in testcaseout json file, read in as list of dictionaries
full_data = {}
base_data = {}
# print(term_path_list)
for termpath in term_path_list:
    infile_name = f"ortho_RPP_inter{termpath[0]}.json"
    infile = open(infile_name)
    base_data = json.load(infile)
    full_data[f"{termpath[0]}"] = ["", ""]
    full_data[f"{termpath[0]}"][0] = termpath[1]
    # print(f"DBUG1: {termpath[0]}: {termpath[1]}")
    full_data[f"{termpath[0]}"][1] = copy.copy(base_data)
    infile.close()
# take in testframeout, read in as pandas file
infile_name = "ortho_DF_inter.csv"
infile = open(infile_name)
dfb = pd.read_csv(infile, index_col=0)
infile.close()

# Next, loop through the dictionaries you want, grabbing the appropriate row by index lookup(?) in the data frame,
# and attaching them to their Uniprot ID references, turning the strings that served as the Data Frame's values
# into lists
# new_dict = {}

dfb.insert(0, "Pathway", "BadPath")
dfb.insert(1, "Pathway_ID", "BadID")
dfb.insert(2, "Reaction", "BadRxn")
dfb.insert(3, "Reaction_ID", "BadID")
dfb.insert(4, "MSU_ID", "")
dfb.insert(5, "RAP_ID", "")

for path in path_list:
    if path in full_data:
        # new_dict[path] = {}
        for rxn in full_data[path][1]:
            # new_dict[path][rxn] = {}
            for prot in full_data[path][1][rxn][1]:
                # new_dict[path][rxn][prot] = {}
                if prot in dfb.index:
                    dfb.loc[prot, "Pathway"] = full_data[path][0]
                    # print(f"DBUG3: {full_data[path][0]}")
                    dfb.loc[prot, "Pathway_ID"] = path
                    # print(f"DBUG4: {path}")
                    dfb.loc[prot, "Reaction"] = full_data[path][1][rxn][0].capitalize()
                    # print(f"DBUG5: {full_data[path][1][rxn][0]}")
                    dfb.loc[prot, "Reaction_ID"] = rxn
                    # print(f"DBUG6: {rxn}")
                    dfb.loc[prot, "MSU_ID"] = full_data[path][1][rxn][1][prot][0]
                    dfb.loc[prot, "RAP_ID"] = full_data[path][1][rxn][1][prot][1]
                print(f"{path} - {rxn} - {prot} processed")
            print(f"{path} - {rxn} processed")
        print(f"{path} processed")

dfb.reset_index(inplace=True)
dfb.rename(columns={'index': 'UniProt_ID'}, inplace=True)
# Finally, use this new dictionary to export into whatever file formats you like.
# df = pd.DataFrame.from_dict(data=dfb,
#                            orient='columns')
dfb.set_index(keys=['Pathway_ID', 'Reaction_ID', 'UniProt_ID'], inplace=True)
dfb.sort_index(axis=0, level=[0, 1, 2], inplace=True)
dfb.to_csv(path_or_buf="ortho_DF_full.csv",
           mode="w")

# dfb.sort_values(by=["Pathway", "Reaction"], inplace=True)

# dfWork = dfb.truncate(before="Actinidia chinensis", axis=1)
df_filtered = dfb.filter(items=species_list)
df_filtered.to_csv(path_or_buf="ortho_DF_filtered.csv", mode="w")

dfWork = df_filtered.applymap(lambda x: len(x.split("|")) if isinstance(x, str) else 0)

# Sorting to pick out which proteins have no orthologs


dfWork.to_csv(path_or_buf="ortho_DF_work.csv", mode="w")
'''if CountFlag is True:
    for i in dfWork.index:
        for c in dfWork.columns:
            if type(dfWork.loc[i, c]) is str:
                tempString = dfWork.loc[i, c]
                tempList = tempString.split("|")
                length = len(tempList)
                dfWork.loc[i, c] = length
            else:
                dfWork.loc[i, c] = 0'''

# Create the various sorted dataframes; use .sum() to find total number of orthologs if ratioFlag is false, or .mean()
# to find the average value, which corresponds to the relative frequency of orthologs compared to rice proteins in a
# given reaction or pathway; this helps even out the data between reactions or pathways with many proteins to begin with
# and those with very few, helping shed light on when a species actually has a relatively high number of orthologs.
if ratioFlag is False:
    df_W = dfWork.groupby(level=0).sum()
    df_R = dfWork.groupby(level=1).sum()
    df_P = dfWork.groupby(level=2).sum()
else:
    df_W = dfWork.groupby(level=0).mean()
    df_R = dfWork.groupby(level=1).mean()
    df_P = dfWork.groupby(level=2).mean()

# Here we find and then filter out those values for which no orthologs are currently present in the database for any
# species.
dfSortW = df_W.apply(np.sum, axis=1)
dfSortW[dfSortW == 0].to_csv(path_or_buf="ortho_DF_sortP.csv", mode="w", header=True)
print("Successfully saved Data Missing Pathway Table")
dfSortR = df_R.apply(np.sum, axis=1)
dfSortR[dfSortR == 0].to_csv(path_or_buf="ortho_DF_sortR.csv", mode="w", header=True)
print("Successfully saved Data Missing Reaction Table")
dfSortP = df_P.apply(np.sum, axis=1)
dfSortP[dfSortP == 0].to_csv(path_or_buf="ortho_DF_sortP.csv", mode="w", header=True)
print("Successfully saved Data Missing Protein Table")

if logFlag is True:
    df_Wm = df_W.drop(labels=dfSortW[dfSortW == 0].index, axis=0,)
    df_Rm = df_R.drop(labels=dfSortR[dfSortR == 0].index, axis=0,)
    df_Pm = df_P.drop(labels=dfSortP[dfSortP == 0].index, axis=0,)

    # Here we take the log of the values; suspect this may be most useful with ratioFlag set to True.
    df_W = np.log2(df_Wm)
    df_W.replace(to_replace=float('-Inf'), value=-1, inplace=True)
    df_R = np.log2(df_Rm)
    df_R.replace(to_replace=float('-Inf'), value=-1, inplace=True)
    df_P = np.log2(df_Pm)
    df_P.replace(to_replace=float('-Inf'), value=-1, inplace=True)
else:
    df_W.drop(labels=dfSortW[dfSortW == 0].index, axis=0, inplace=True)
    df_R.drop(labels=dfSortR[dfSortR == 0].index, axis=0, inplace=True)
    df_P.drop(labels=dfSortP[dfSortP == 0].index, axis=0, inplace=True)
kwargs = {"family": 'arial', "style": 'italic'}
# Pathway heatmap and clustermap
if levelFlag is 'A' or levelFlag is 'W':
    aWx = sns.heatmap(data=df_W, cmap="viridis", xticklabels=True,
                      # yticklabels=True,
                      )
    aWx.set_xticklabels(aWx.get_xticklabels(), fontsize=5, **kwargs)
    plt.savefig(f"heatMap_TypeWT.pdf", bbox_inches='tight')
    print("Successfully saved Pathway Heatmap")

    bWx = sns.clustermap(data=df_W, cmap="viridis",
                         xticklabels=True,
                         yticklabels=True,
                         col_cluster=False,
                         figsize=[10, 24],
                         cbar_kws={'label': 'Gene-Orthology-Based Pathway Projections',
                                   'orientation': 'horizontal'},
                         )
    bWx.ax_col_dendrogram.set_visible(False)
    W_DendroBox = bWx.ax_col_dendrogram.get_position()
    print(W_DendroBox)
    W_DendroBox.y0 =(W_DendroBox.y0 + 9 * W_DendroBox.y1) / 10
    W_DendroWid = W_DendroBox.y0 - W_DendroBox.y1
    W_DendroBox.y0 = bWx.ax_col_dendrogram.get_position().y0
    W_DendroBox.y1 = W_DendroBox.y0 - W_DendroWid
    print(W_DendroBox)
    bWx.cax.set_position(W_DendroBox)
    bWx.cax.xaxis.set_ticks_position(position="top")
    bWx.cax.xaxis.set_label_position(position="top")
    bWx.ax_heatmap.tick_params(axis='x', labelsize=7)
    bWx.ax_heatmap.tick_params(axis='y', labelsize=7)
    bWx.ax_heatmap.set_ylabel("Pathway stIDs")
    bWx.ax_heatmap.set_xticklabels(bWx.ax_heatmap.get_xticklabels(), **kwargs)
    plt.savefig(f"clusterMap_TypeWT.pdf", bbox_inches='tight')
    print("Successfully saved Pathway Clustermap")
    plt.show()

# Reaction heatmap and clustermap
if levelFlag is 'A' or levelFlag is 'R':
    aRx = sns.heatmap(data=df_R, cmap="viridis", xticklabels=True,
                      # yticklabels=True,
                      )
    aRx.set_xticklabels(aRx.get_xmajorticklabels(), fontsize=4, **kwargs)
    plt.savefig(f"heatMap_TypeRT.png", bbox_inches='tight')
    print("Successfully saved Reaction Heatmap")

    bRx = sns.clustermap(data=df_R, cmap="viridis",
                         xticklabels=True,
                         yticklabels=2,
                         col_cluster=False,
                         figsize=[10, 24],
                         cbar_kws={'label': 'Gene-Orthology-Based Reaction Projections',
                                   'orientation': 'horizontal'},
                         )
    bRx.ax_col_dendrogram.set_visible(False)
    R_DendroBox = bRx.ax_col_dendrogram.get_position()
    print(R_DendroBox)
    R_DendroBox.y0 = (R_DendroBox.y0 + 9 * R_DendroBox.y1) / 10
    R_DendroWid = R_DendroBox.y0 - R_DendroBox.y1
    R_DendroBox.y0 = bRx.ax_col_dendrogram.get_position().y0
    R_DendroBox.y1 = R_DendroBox.y0 - R_DendroWid
    print(R_DendroBox)
    bRx.cax.set_position(R_DendroBox)
    bRx.cax.xaxis.set_ticks_position(position="top")
    bRx.cax.xaxis.set_label_position(position="top")
    bRx.ax_heatmap.tick_params(axis='x', labelsize=6)
    bRx.ax_heatmap.tick_params(axis='y', labelsize=5)
    bRx.ax_heatmap.set_ylabel("Reaction stIDs")
    bRx.ax_heatmap.set_xticklabels(bRx.ax_heatmap.get_xticklabels(), **kwargs)
    plt.savefig(f"clusterMap_TypeRT.png", bbox_inches='tight')
    print("Successfully saved Reaction Clustermap")
    plt.show()

# Protein heatmap and clustermap
if levelFlag is 'A' or levelFlag is 'P':
    aPx = sns.heatmap(data=df_P, cmap="viridis",
                      # xticklabels=True,
                      # yticklabels=True,
                      )
    aPx.set_xticklabels(aPx.get_xmajorticklabels(), fontsize=5, **kwargs)
    plt.savefig(f"heatMap_TypePT.png", bbox_inches='tight')
    print("Successfully saved Protein Heatmap")

    bPx = sns.clustermap(data=df_P, cmap="viridis",
                         xticklabels=True,
                         yticklabels=6,
                         col_cluster=False,
                         figsize=[10, 24],
                         cbar_kws={'label': 'Gene-Orthology-Based Protein Projections',
                                   'orientation': 'horizontal'},
                         )
    bPx.ax_col_dendrogram.set_visible(False)
    P_DendroBox = bPx.ax_col_dendrogram.get_position()
    print(P_DendroBox)
    P_DendroBox.y0 = (P_DendroBox.y0 + 9 * P_DendroBox.y1) / 10
    P_DendroWid = P_DendroBox.y0 - P_DendroBox.y1
    P_DendroBox.y0 = bPx.ax_col_dendrogram.get_position().y0
    P_DendroBox.y1 = P_DendroBox.y0 - P_DendroWid
    print(P_DendroBox)
    bPx.cax.set_position(P_DendroBox)
    bPx.cax.xaxis.set_ticks_position(position="top")
    bPx.cax.xaxis.set_label_position(position="top")
    bPx.ax_heatmap.tick_params(axis='x', labelsize=6)
    bPx.ax_heatmap.tick_params(axis='y', labelsize=5)
    bPx.ax_heatmap.set_ylabel("UniProt IDs in Oryza Sativa")
    bPx.ax_heatmap.set_xticklabels(bPx.ax_heatmap.get_xticklabels(), **kwargs)
    plt.savefig(f"clusterMap_TypePT.png", bbox_inches='tight')
    print("Successfully saved Protein Clustermap")
    plt.show()
print("Successfully reached end of program!")
