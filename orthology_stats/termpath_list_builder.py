import json
import requests
import sys
import pandas as pd

url_base = 'https://plantreactomedev.gramene.org/ContentService'
headers = {'accept': 'application/json'}


def recurs_get_paths(sub_dict, termDict, path_list):
    pathTypes = ['Pathway', 'TopLevelPathway']
    rxnTypes = ['Reaction', 'BlackBoxEvent']
    if sub_dict['stId'] in path_list and sub_dict['type'] in pathTypes:
        term_path_flag = False
        for child in sub_dict['children']:
            if child['type'] == 'Pathway':
                path_list.append(child['stId'])
                recurs_get_paths(child, termDict, path_list)
            elif child['type'] in rxnTypes and term_path_flag is False:
                term_path_flag = True
                if f"{sub_dict['stId']}" not in termDict:
                    termDict[f"{sub_dict['stId']}"] = False
    elif sub_dict['type'] in pathTypes:
        for child in sub_dict['children']:
            if child['type'] in pathTypes:
                recurs_get_paths(child, termDict, path_list)
    return


if __name__ == '__main__':
    sys.argv = ['termpath_list_builder.py',
                'R-OSA-2744345']
    if len(sys.argv) > 1:
        entryId = sys.argv[1:]
    else:
        entryId = ['R-OSA-2744345']
    path_url = '/data/eventsHierarchy/4530' # rice as reference (NCBI tax id)
    full_url = url_base + path_url
    response = requests.get(full_url, headers=headers).json()
    base_dict = {}
    term_dict = {'Header': False}
    for TopLevel in response:
        if TopLevel['name'] == "Metabolism and regulation":  # configure top level pathway by name here
            base_dict = TopLevel
            break
    recurs_get_paths(base_dict, term_dict, entryId)
    df = pd.DataFrame.from_dict(data=term_dict, orient='index')

    df.to_csv(path_or_buf="termpath_checklist.csv", mode="w", header="Completed?")
    # with open(term, 'rb+') as filehandle:
    #    filehandle.seek(-1, os.SEEK_END)
    #    filehandle.truncate()
    # outfile = open(outDict, "a")
    # outfile.write("\n}")
    # outfile.close()
