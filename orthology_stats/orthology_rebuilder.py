#!C:\Program Files (x86)\Python
import json
from typing import Any, Union
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame
from pandas.io.parsers import TextFileReader

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

# if __name__ == '__main__':
CountFlag = True

# something something interface for user to decide which pathways they want?
path_lists = ["R-OSA-1119263"]
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
dfb.insert(1, "Reaction", "BadRxn")
dfb.insert(2, "MSU_ID", "")
dfb.insert(3, "RAP_ID", "")

for path in path_lists:
    if path in base_data:
        # new_dict[path] = {}
        for rxn in base_data[path]:
            # new_dict[path][rxn] = {}
            for prot in base_data[path][rxn]:
                # new_dict[path][rxn][prot] = {}
                if prot in dfb.index:
                    dfb.loc[prot, "Pathway"] = path
                    dfb.loc[prot, "Reaction"] = rxn.capitalize()
                    if "MSU_ID" in base_data[path][rxn][prot]:
                        dfb.loc[prot,"MSU_ID"] = base_data[path][rxn][prot]["MSU_ID"]
                    if "RAP_ID" in base_data[path][rxn][prot]:
                        dfb.loc[prot, "RAP_ID"] = base_data[path][rxn][prot]["RAP_ID"]

dfb.sort_values(by=["Pathway", "Reaction"],
                inplace=True)

# Finally, use this new dictionary to export into whatever file formats you like.
df = pd.DataFrame.from_dict(data=dfb,
                            orient='columns')
df.to_csv(path_or_buf="ortho_DF_full.csv",
          mode="w")

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


ax = sns.heatmap(data=dfWork,
                 cmap="BuGn")
plt.savefig("heatMapExample.png", bbox_inches='tight')

bx = sns.clustermap(data=dfWork,
                    cmap="BuGn")
plt.savefig("clusterMapExample.png", bbox_inches='tight')
