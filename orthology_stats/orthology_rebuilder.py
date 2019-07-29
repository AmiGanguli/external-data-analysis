#!C:\Program Files (x86)\Python
import json
import requests
import sys
from typing import Any, Union
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame
from pandas.io.parsers import TextFileReader

url_base = 'https://plantreactomedev.gramene.org/ContentService'
headers = {'accept': 'application/json'}

species_list=["Actinidia chinensis", "Aegilops tauschii", "Amborella trichopoda", "Arabidopsis halleri",
              "Arabidopsis lyrata", "Arabidopsis thaliana", "Arachis duranensis", "Arachis ipaensis", "Beta vulgaris",
              "Brachypodium distachyon", "Brassica napus", "Brassica oleracea", "Brassica rapa", "Cajanus cajan",
              "Capsicum annuum", "Chlamydomonas reinhardtii", "Chondrus crispus", "Cicer arietinum", "Citrus sinensis",
              "Coffea canephora", "Corchorus capsularis", "Cucumis sativus", "Cyanidioschyzon merolae", "Daucus carota",
              "Dioscorea rotundata", "Erythranthe guttata", "Eucalyptus grandis", "Fragaria vesca",
              "Galdieria sulphuraria", "Glycine max", "Gossypium raimondii", "Helianthus annuus", "Hordeum vulgare",
              "Jatropha curcas", "Leersia perrieri", "Lupinus angustifolius", "Malus domestica", "Manihot esculenta",
              "Medicago truncatula", "Musa acuminata", "Nicotiana attenuata", "Oryza australiensis", "Oryza barthii",
              "Oryza brachyantha", "Oryza glaberrima", "Oryza glumaepatula", "Oryza longistaminata",
              "Oryza meridionalis", "Oryza meyeriana var. granulata", "Oryza minuta", "Oryza nivara",
              "Oryza officinalis", "Oryza punctata", "Oryza rufipogon", "Oryza sativa Indica Group",
              "Oryza sativa aus subgroup", "Oryza sativa subsp. japonica", "Ostreococcus lucimarinus",
              "Panicum hallii FIL2", "Panicum hallii var. hallii HAL2", "Phaseolus vulgaris", "Phoenix dactylifera",
              "Physcomitrella patens", "Picea abies", "Pinus taeda", "Populus trichocarpa", "Prunus persica",
              "Selaginella moellendorffii", "Setaria italica", "Solanum lycopersicum", "Solanum tuberosum",
              "Sorghum bicolor", "Synechocystis sp. PCC 6803", "Theobroma cacao", "Trifolium pratense",
              "Triticum aestivum", "Triticum dicoccoides", "Triticum turgidum", "Triticum urartu", "Vigna angularis",
              "Vigna radiata", "Vitis vinifera", "Zea mays"]

def recurs_get_paths(sub_dict, path_list, IDs):
    pathTypes = ['Pathway', 'TopLevelPathway']
    rxnTypes = ['Reaction', 'BlackBoxEvent']
    if sub_dict['stId'] in path_list and sub_dict['type'] in pathTypes:
        for child in sub_dict['children']:
            if child['type'] == 'Pathway':
                path_list.append(child['stId'])
                IDs[f"{child['stId']}"] = child['displayName']
                recurs_get_paths(child, path_list, IDs)
            elif child['type'] in rxnTypes:
                IDs[f"{child['displayName']}"] = child['stId']
    elif sub_dict['type'] in pathTypes:
        for child in sub_dict['children']:
            if child['type'] in pathTypes:
                recurs_get_paths(child, path_list, IDs)
    return


# if __name__ == '__main__':
CountFlag = True

# something something interface for user to decide which pathways they want?
sys.argv = ['orthology_rebuilder.py', 'P',
            'R-OSA-1119263']
# P -> Protein Level; R -> Reaction Level; W -> pathWay Level
levelFlag = sys.argv[1]
path_list = sys.argv[2:]

path_url = '/data/eventsHierarchy/4530'
full_url = url_base+path_url
response = requests.get(full_url, headers=headers).json()

base_dict = {}
for TopLevel in response:
    if TopLevel['name'] == "Metabolism and regulation":
        base_dict = TopLevel
        break
ID_dict = {}
recurs_get_paths(base_dict, path_list, ID_dict)


# take in testcaseout json file, read in as list of dictionaries
infile_name = "ortho_RPP_inter.json"
infile = open(infile_name)
base_data = json.load(infile)
# take in testframeout, read in as pandas file
infile_name = "ortho_DF_inter.csv"
infile = open(infile_name)
dfb = pd.read_csv(infile, index_col=0)

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
        for rxn in base_data[path]:
            # new_dict[path][rxn] = {}
            for prot in base_data[path][rxn]:
                # new_dict[path][rxn][prot] = {}
                if prot in dfb.index:
                    dfb.loc[prot, "Pathway"] = ID_dict[path]
                    dfb.loc[prot, "Pathway_ID"] = path
                    dfb.loc[prot, "Reaction"] = rxn.capitalize()
                    dfb.loc[prot, "Reaction_ID"] = ID_dict[rxn]
                    if "MSU_ID" in base_data[path][rxn][prot]:
                        dfb.loc[prot,"MSU_ID"] = base_data[path][rxn][prot]["MSU_ID"]
                    if "RAP_ID" in base_data[path][rxn][prot]:
                        dfb.loc[prot, "RAP_ID"] = base_data[path][rxn][prot]["RAP_ID"]

dfb.reset_index(inplace=True)
dfb.rename_axis(columns={'index':'UniProt_ID'}, inplace=True)
dfb.set_index(keys=['Pathway', 'Reaction', 'UniProt_ID'],inplace=True)


# Finally, use this new dictionary to export into whatever file formats you like.
df = pd.DataFrame.from_dict(data=dfb,
                            orient='columns')
df.to_csv(path_or_buf="ortho_DF_full.csv",
          mode="w")

# dfb.sort_values(by=["Pathway", "Reaction"], inplace=True)
if levelFlag is 'W':
    dfb.groupby(level=0)
elif levelFlag is 'R':
    dfb.groupby(level=1)
elif levelFlag is 'P':
    dfb.groupby(level=2)
else:
    print("Sorry, something's gone wrong; you may have entered what you wanted to sort by incorrectly.\nCurrently "
          "acceptable choices are:\n\'P\' - UniProtID\n\'R\' - Reaction\n\'W\' - Terminal Pathway")
    exit()
# dfWork = dfb.truncate(before="Actinidia chinensis", axis=1)
dfWork = dfb.filter(items=species_list)

if CountFlag is True:
    for i in dfWork.index:
        for c in dfWork.columns:
            if type(dfWork.loc[i, c]) is str:
                tempString = dfWork.loc[i, c]
                tempList = tempString.split("|")
                length = len(tempList)
                dfWork.loc[i, c] = length
            else:
                dfWork.loc[i, c] = 0


ax = sns.heatmap(data=dfWork, cmap="BuGn", figsize=(10, 16))
plt.savefig(f"heatMap_Type{levelFlag}.png", bbox_inches='tight')

bx = sns.clustermap(data=dfWork, cmap="BuGn", figsize=(16, 10))
plt.savefig(f"clusterMap_Type{levelFlag}.png", bbox_inches='tight')
