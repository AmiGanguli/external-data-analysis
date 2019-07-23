#!C:\Program Files (x86)\Python
import json
from typing import Any, Union

import pandas as pd
from pandas import DataFrame
from pandas.io.parsers import TextFileReader

#if __name__ == '__main__':


# something something interface for user to decide which pathways they want?
path_lists = []
# take in testcaseout json file, read in as list of dictionaries
infile = "testcaseout.txt"
base_data = json.load(infile)
# take in testframeout, read in as pandas file
infile_2 = "testframeout.txt"
dfb = pd.read_csv(infile_2, index_col=0)

# Next, loop through the dictionaries you want, grabbing the appropriate row by index lookup(?) in the data frame,
# and attaching them to their Uniprot ID references, turning the strings that served as the Data Frame's values
# into lists
# new_dict = {}

for path in path_lists:
    if path in base_data:
        # new_dict[path] = {}
        for rxn in base_data[path]:
            # new_dict[path][rxn] = {}
            for prot in base_data[path]:
                # new_dict[path][rxn][prot] = {}
                if prot in dfb.index:
                    curRow = dfb.loc[prot]
                    for species in dfb.columns:
                        orthCur = curRow.loc[species]
                        # new_dict[path][rxn][prot][species] = orthCur.split(',')

# Finally, use this new dictionary to export into whatever file formats you like.
print(dfb)