#!C:\Program Files (x86)\Python
import json
import sys
import fastcluster
import matplotlib.pyplot as plt
import pandas as pd
import requests
import seaborn as sns
import copy

url_base = 'https://plantreactomedev.gramene.org/ContentService'
headers = {'accept': 'application/json'}

species_list = ["Actinidia chinensis", "Aegilops tauschii", "Amborella trichopoda", "Arabidopsis halleri",
                "Arabidopsis lyrata", "Arabidopsis thaliana", "Arachis duranensis", "Arachis ipaensis", "Beta vulgaris",
                "Brachypodium distachyon", "Brassica napus", "Brassica oleracea", "Brassica rapa", "Cajanus cajan",
                "Capsicum annuum", "Chlamydomonas reinhardtii", "Chondrus crispus", "Cicer arietinum",
                "Citrus sinensis", "Coffea canephora", "Corchorus capsularis", "Cucumis sativus",
                "Cyanidioschyzon merolae", "Daucus carota", "Dioscorea rotundata", "Erythranthe guttata",
                "Eucalyptus grandis", "Fragaria vesca", "Galdieria sulphuraria", "Glycine max", "Gossypium raimondii",
                "Helianthus annuus", "Hordeum vulgare", "Jatropha curcas", "Leersia perrieri", "Lupinus angustifolius",
                "Malus domestica", "Manihot esculenta", "Medicago truncatula", "Musa acuminata", "Nicotiana attenuata",
                "Oryza australiensis", "Oryza barthii", "Oryza brachyantha", "Oryza glaberrima", "Oryza glumaepatula",
                "Oryza longistaminata", "Oryza meridionalis", "Oryza meyeriana var. granulata", "Oryza minuta",
                "Oryza nivara", "Oryza officinalis", "Oryza punctata", "Oryza rufipogon", "Oryza sativa Indica Group",
                "Oryza sativa aus subgroup", "Oryza sativa subsp. japonica", "Ostreococcus lucimarinus",
                "Panicum hallii FIL2", "Panicum hallii var. hallii HAL2", "Phaseolus vulgaris", "Phoenix dactylifera",
                "Physcomitrella patens", "Picea abies", "Pinus taeda", "Populus trichocarpa", "Prunus persica",
                "Selaginella moellendorffii", "Setaria italica", "Solanum lycopersicum", "Solanum tuberosum",
                "Sorghum bicolor", "Synechocystis sp. PCC 6803", "Theobroma cacao", "Trifolium pratense",
                "Triticum aestivum", "Triticum dicoccoides", "Triticum turgidum", "Triticum urartu", "Vigna angularis",
                "Vigna radiata", "Vitis vinifera", "Zea mays"]


def recurs_get_paths(sub_dict, path_list, term_path_list, IDs):
    pathTypes = ['Pathway', 'TopLevelPathway']
    rxnTypes = ['Reaction', 'BlackBoxEvent']
    if sub_dict['stId'] in path_list and sub_dict['type'] in pathTypes:
        # termFlag = False
        for child in sub_dict['children']:
            if child['type'] == 'Pathway':
                path_list.append(child['stId'])
                IDs[f"{child['stId']}"] = child['name']

                recurs_get_paths(child, path_list, IDs)
            elif child['type'] in rxnTypes:
                # termFlag = True
                sub_term_list = ["", ""]
                sub_term_list[0] = sub_dict['stId']
                sub_term_list[1] = sub_dict['name']
                term_path_list.append(sub_term_list)
                IDs[f"{child['name']}"] = child['stId']
    elif sub_dict['type'] in pathTypes:
        for child in sub_dict['children']:
            if child['type'] in pathTypes:
                recurs_get_paths(child, path_list, IDs)
    return


# if __name__ == '__main__':
CountFlag = True

# something something interface for user to decide which pathways they want?
sys.argv = ['orthology_rebuilder.py', 'P', '-r',
            'R-OSA-2744345']
# P -> Protein Level; R -> Reaction Level; W -> pathWay Level
levelFlag = sys.argv[1]
path_list = sys.argv[3:]
term_path_list = []
ratioFlag = False
if sys.argv[2] is '-r':
    ratioFlag = True

path_url = '/data/eventsHierarchy/4530'
full_url = url_base + path_url
response = requests.get(full_url, headers=headers).json()

base_dict = {}
for TopLevel in response:
    if TopLevel['name'] == "Metabolism and regulation":
        base_dict = TopLevel
        break
ID_dict = {}
recurs_get_paths(base_dict, path_list, term_path_list, ID_dict)

# take in testcaseout json file, read in as list of dictionaries
full_data = {}
base_data = {}
for termpath in term_path_list:
    infile_name = f"ortho_RPP_inter{termpath}.json"
    infile = open(infile_name)
    base_data = json.load(infile)
    full_data[f"{termpath[0]}"] = ["", ""]
    full_data[f"{termpath[0]}"][0] = termpath[1]
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
    if path in base_data:
        # new_dict[path] = {}
        for rxn in base_data[path][1]:
            # new_dict[path][rxn] = {}
            for prot in base_data[path][1][rxn][1]:
                # new_dict[path][rxn][prot] = {}
                if prot in dfb.index:
                    dfb.loc[prot, "Pathway"] = base_data[path][0]
                    dfb.loc[prot, "Pathway_ID"] = path
                    dfb.loc[prot, "Reaction"] = base_data[path][1][rxn][0].capitalize()
                    dfb.loc[prot, "Reaction_ID"] = rxn
                    dfb.loc[prot, "MSU_ID"] = base_data[path][1][rxn][1][prot][0]
                    dfb.loc[prot, "RAP_ID"] = base_data[path][1][rxn][1][prot][1]

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

if ratioFlag is False:
    if levelFlag is 'W':
        df_med = dfWork.groupby(level=0).sum()
    elif levelFlag is 'R':
        df_med = dfWork.groupby(level=1).sum()
    elif levelFlag is 'P':
        df_med = dfWork.groupby(level=2).sum()
    elif levelFlag is 'A':
        df_W = dfWork.groupby(level=0).sum()
        df_R = dfWork.groupby(level=1).sum()
        df_P = dfWork.groupby(level=2).sum()
    else:
        print("Sorry, something's gone wrong; you may have entered what you wanted to sort by incorrectly.\nCurrently "
              "acceptable choices are:\n\'P\' - UniProtID\n\'R\' - Reaction\n\'W\' - Terminal Pathway")
        exit()
else:
    if levelFlag is 'W':
        df_med = dfWork.groupby(level=0).mean()
    elif levelFlag is 'R':
        df_med = dfWork.groupby(level=1).mean()
    elif levelFlag is 'P':
        df_med = dfWork.groupby(level=2).mean()
    elif levelFlag is 'A':
        df_W = dfWork.groupby(level=0).mean()
        df_R = dfWork.groupby(level=1).mean()
        df_P = dfWork.groupby(level=2).mean()
    else:
        print("Sorry, something's gone wrong; you may have entered what you wanted to sort by incorrectly.\nCurrently "
              "acceptable choices are:\n\'P\' - UniProtID\n\'R\' - Reaction\n\'W\' - Terminal Pathway")
        exit()

if levelFlag is 'A':
    aWx = sns.heatmap(data=df_W, cmap="plasma", xticklabels=True, yticklabels=True)
    plt.savefig(f"heatMap_TypeW.png", bbox_inches='tight')

    bWx = sns.clustermap(data=df_W, cmap="plasma", xticklabels=True, yticklabels=True)
    bWx.ax_heatmap.tick_params(axis='x', labelsize=6)
    plt.savefig(f"clusterMap_TypeW.png", bbox_inches='tight')

    aRx = sns.heatmap(data=df_R, cmap="plasma", xticklabels=True, yticklabels=True)
    plt.savefig(f"heatMap_TypeR.png", bbox_inches='tight')

    bRx = sns.clustermap(data=df_R, cmap="plasma", xticklabels=True, yticklabels=True)
    bRx.ax_heatmap.tick_params(axis='x', labelsize=6)
    plt.savefig(f"clusterMap_TypeR.png", bbox_inches='tight')

    aPx = sns.heatmap(data=df_P, cmap="plasma", xticklabels=True, yticklabels=True)
    plt.savefig(f"heatMap_TypeP.png", bbox_inches='tight')

    bPx = sns.clustermap(data=df_P, cmap="plasma", xticklabels=True, yticklabels=True)
    bPx.ax_heatmap.tick_params(labelsize=5)
    plt.savefig(f"clusterMap_TypeP.png", bbox_inches='tight')
else:
    df_med.to_csv(path_or_buf="ortho_DF_med.csv", mode="w")
    ax = sns.heatmap(data=df_med, cmap="plasma")
    plt.savefig(f"heatMap_Type{levelFlag}.png", bbox_inches='tight')

    bx = sns.clustermap(data=df_med, cmap="plasma")
    plt.savefig(f"clusterMap_Type{levelFlag}.png", bbox_inches='tight')
