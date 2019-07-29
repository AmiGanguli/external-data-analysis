#!C:\Program Files (x86)\Python
import requests
import re
import sys
import math
import time
import json
import pandas as pd
import copy
import os
url_base = 'https://plantreactomedev.gramene.org/ContentService'
headers = {'accept': 'application/json'}

#########################################################################
#                     Method Call Hierarchy                             #
#    1) Call Event Hierarchy method on species of interest              #
#           If only one toplevel path is desired, grab here.            #
#    2) get_path_data to recursively search down pathways until         #
#           reaction-like events are found.                             #
#    3) Call Participants on reaction-like events                       #
#           Grab anything that might have EWAS, mostly Cat. Activity    #
#    4) Call Query and search through hasMember                         #
#           If SimpleEntity, ignore                                     #
#           If not an EWAS, call Query recursively until EWAS           #
#    5) Call Query on EWAS and search through inferredTo                #
#           Add orthologs to species dict and reaction dict             #
#########################################################################


def usefulPrint(rxn_dict, path_depth, depth, outfile):
    if depth == 0:
        path_reach = path_depth*"Pathway\t"
        outfile.write(f"{path_reach}Reaction\tUniProt_ID\tMSU?\tRAP?\tSpecies\tOrthologs\n")
    if 'species' in rxn_dict:
        extra = 0
        if 'MSU' in rxn_dict:
            outfile.write(f"{rxn_dict['MSU']}\t")
            extra += 4
        else:
            outfile.write(f"No MSU\t")
            extra += 4
        if 'RAP' in rxn_dict:
            outfile.write(f"{rxn_dict['RAP']}\t")
            extra += 4
        else:
            outfile.write(f"No RAP\t")
            extra +=4
        usefulPrint(rxn_dict['species'], path_depth-1, depth + extra, outfile)
        # outfile.write(f'\n')
        if len(rxn_dict['species']) > 0 and len(rxn_dict) > 0:
            print("no linebreak")
            # for tabs in range(0, depth):
            #    outfile.write(f"\t")
        else:
            outfile.write(f'\n')
    else:
        # rxn_keys = [*rxn_dict]
        for element in rxn_dict:
            print(f"{element} Depth: {depth}")
            outfile.write(f"{element}\t")
            extra = math.floor(len(element)/4) + 1
            print(f"{element} Extra Depth added: {extra}")
            if type(rxn_dict) == dict:
                usefulPrint(rxn_dict[element], path_depth-1, depth + extra, outfile)
            else:
                outfile.write(f'\n')


def get_multi_product_data(setId, df_dict, ortho_spec, UniProtId):
    new_path = f'/data/query/{setId}'
    new_full = url_base + new_path
    new_response = requests.get(new_full, headers=headers).json()
    for member in new_response['hasMember']:
        print(f'DBUG9.2 --- Name: {member["name"][0]} --- Schema Class: {member["schemaClass"]} --- Query: '
              f'Query[hasMember]')
        if member['schemaClass'] == "EntityWithAccessionedSequence":
            if ortho_spec in df_dict:
                if UniProtId not in df_dict[ortho_spec]:
                    df_dict[ortho_spec][UniProtId] = []
                df_dict[ortho_spec][UniProtId].append(member['name'][0])
    '''if new_response.schemaClass == "DefinedSet":
        for member in new_response.hasMember:
            if member.schemaClass != "SimpleEntity":
                get_multi_product_data(member.stId,rxn_dict)
    elif new_response.schemaClass == "Complex":
        for component in new_response.hasComponent:
            if component.schemaClass != "SimpleEntity:
                get_multi_product_data(member.stId, rxn_dict)
            elif new_response.schemaClass == "EntityWithAccessionSequence":
                get_multi_ortho_data(new_response, rxn_dict)'''
    return


def get_ortho_data(ewas, df_dict, UniProtId):
    if 'inferredTo' in ewas:
        for ortholog in ewas['inferredTo']:
            ortho_spec = ortholog['speciesName']
            print(f'DBUG9.1 --- Name: {ortholog["name"][0]} --- Schema Class: {ortholog["schemaClass"]} --- Query: '
                  f'Query[inferredTo]')

            if ortholog['schemaClass'] == "EntityWithAccessionedSequence":
                if ortho_spec in df_dict:
                    df_dict[ortho_spec][UniProtId].append(ortholog['name'][0])
            elif ortholog['schemaClass'] == 'DefinedSet':
                get_multi_product_data(ortholog['stId'], df_dict, ortho_spec, UniProtId)
    return


def get_prot_data(ewas, rxn_dict, df_dict):
    if 'identifier' in ewas['referenceEntity']:
        UniProtId = ewas['referenceEntity']['identifier']
    else:
        UniProtId = ewas['referenceEntity']['secondaryIdentifier'][0]

    rxn_dict[UniProtId] = {}
    print(f'DBUG8.1 --- Name: {ewas["displayName"]} --- stId: {ewas["stId"]} --- UniProt: {UniProtId} --- Query: Query')

    if 'geneName' in ewas['referenceEntity']:
        RAP_flag = False
        MSU_flag = False
        for name in ewas['referenceEntity']['geneName']:
            match = re.match("OS..G........", name.upper())
            if match and RAP_flag is False:
                rxn_dict[UniProtId]["RAP_ID"] = name[0:12]
                print(f'DBUG8.2 --- RAP ID found: {name[0:12]}')
                RAP_flag = True
            match1 = re.match("LOC_OS..G.....", name.upper())
            if match1 and MSU_flag is False:
                rxn_dict[UniProtId]["MSU_ID"] = name
                print(f'DBUG8.3 --- MSU ID found: {name}')
                MSU_flag = True
            if RAP_flag is True and MSU_flag is True:
                print(f'DBUG8.4 --- both Secondary IDs found')
                break

    for species in df_dict:
        if UniProtId not in df_dict[species]:
            df_dict[species][UniProtId] = []
    get_ortho_data(ewas, df_dict, UniProtId)

    for species in df_dict:
        if len(df_dict[species][UniProtId]) == 0:
            df_dict[species][UniProtId].append("")
        temp_string = "|".join(df_dict[species][UniProtId])
        df_dict[species][UniProtId] = copy.copy(temp_string)


def get_product_data(entityId, rxn_dict, df_dict):
    new_path = f'/data/query/{entityId}'
    new_full = url_base + new_path
    new_response = requests.get(new_full, headers=headers).json()

    if 'stId' in new_response:
        newId = new_response['stId']
    else:
        newId = new_response['dbId']

    print(f'DBUG6.1 --- Name: {new_response["displayName"]} --- ID: {newId} --- Schema Class: '
          f'{new_response["schemaClass"]} --- Query: Query')

    if new_response['schemaClass'] == "CatalystActivity":
        get_product_data(new_response['physicalEntity']['stId'], rxn_dict, df_dict)

    elif new_response['schemaClass'] == "DefinedSet":
        for member in new_response['hasMember']:
            if member['schemaClass'] != "SimpleEntity":
                print(f'DBUG7.2 --- Name: {member["displayName"]} --- stID: {member["stId"]} --- Schema Class: '
                      f'{member["schemaClass"]} --- Query: Query')
                get_product_data(member['stId'], rxn_dict, df_dict)

    elif new_response['schemaClass'] == 'Complex':
        for component in new_response['hasComponent']:
            if component['schemaClass'] != "SimpleEntity":
                print(f'DBUG7.3 --- Name: {component["displayName"]} --- stID: {component["stId"]} --- Schema Class: '
                      f'{component["schemaClass"]} --- Query: Query')
                get_product_data(component['stId'], rxn_dict, df_dict)

    elif new_response['schemaClass'] == "EntityWithAccessionedSequence":
        get_prot_data(new_response, rxn_dict, df_dict)


# altered version of get_parts_data that skips over CatalystActivity call
def get_parts_data(event, rxn_dict, df_dict):
    new_path = f'/data/participants/{event["stId"]}/participatingPhysicalEntities'
    new_full = url_base + new_path
    new_response = requests.get(new_full, headers=headers).json()

    for party in new_response:
        sClass = party['schemaClass']

        if sClass == "CatalystActivity":
            print(f'DBUG5.1 --- Name: {party["displayName"]} --- dbId: {party["peDbId"]} --- Query: '
                  f'participatingPhysicalEntities')
            get_product_data(party['peDbId'], rxn_dict, df_dict)

        elif sClass == "DefinedSet" or sClass == "Complex" or sClass == "EntityWithAccessionedSequence":
            print(f'DBUG5.2 --- Name: {party["displayName"]} --- stId: {party["stId"]} --- Query: '
                  f'participatingPhysicalEntities')
            get_product_data(party['stId'], rxn_dict, df_dict)


# Something of a wrapper function to reach CatalystActivity
'''def get_parts_data(event, rxn_dict):
    new_path = f'/data/participants/{event["stId"]})'
    new_full = url_base + new_path
    new_response = requests.get(new_full, headers=headers).json()
    for party in new_response:
        if party['schemaClass'] == "CatalystActivity":
            print(f'DBUG6: {party}')
            get_product_data(party['peDbId'], rxn_dict)
    return'''


def build_species_dict(df_dict):
    new_path = f'/data/species/all'
    new_full = url_base + new_path
    spec_res = requests.get(new_full, headers=headers).json()

    for species in spec_res:
        if species["displayName"] not in df_dict:
            df_dict[species["displayName"]] = {}


def term_path_adapter(sub_dict):
    print(f'DBUG4.1 --- Terminal Path Reached --- Name: {sub_dict["name"]} --- stID: {sub_dict["stId"]}')

    rxn_dict = {}
    df_dict = {}
    build_species_dict(df_dict)

    for child in sub_dict['children']:
        print(f'DBUG4.2 --- Name: {child["name"]} --- stID: {child["stId"]} --- Type: {child["type"]} --- '
              f'Query: eventsHierarchy subtree')
        rxn_dict[child['name']] = {}
        get_parts_data(child, rxn_dict[child['name']], df_dict)

    start_time1 = time.time()

    outfile = open(outDict, "a")
    entry = (f"\n\"{sub_dict['stId']}\": " + json.dumps(rxn_dict, indent=4) + ",")
    outfile.write(entry)
    outfile.close()
    time2 = time.time() - start_time1
    print("--- %s seconds ---" % time2)

    df = pd.DataFrame.from_dict(data=df_dict, orient='columns')
    df.to_csv(path_or_buf=outFrame, mode="a", )
    df_dict.clear()
    rxn_dict.clear() # Not sure if this is the best way to free the dictionary's memory


# TODO : simplify tree walk
# Dives into pathways and constructs pathway hierarchy
def get_path_data(sub_dict, path_list):
    rxn_flag = False
    for child in sub_dict['children']:
        # print(f'3.5: {child}')
        if child['stId'] in path_list:
            if 'children' in child:
                for grandkid in child['children']:
                    if grandkid['type'] == 'Pathway':
                        path_list.append(grandkid['stId'])
                        print(f"DBUG_CHILD: Child reached, grandchildren ID grabbed: {path_list}\t---\t{grandkid}")
        if child['stId'] not in path_list:
            if 'children' in child:
                grandflag = False
                for grandkid in child['children']:
                    if grandkid['stId'] in path_list:
                        if child['stId'] not in path_list:
                            path_list.append(child['stId'])
                            print(f"DBUG_GRAND: Grandchild desired, child ID grabbed\t---\t{path_list}")
                            grandflag = True
                            break
                if grandflag is False:
                    for grandkid in child['children']:
                        if 'children' in grandkid:
                            for greatgrandkid in grandkid['children']:
                                if greatgrandkid['stId'] in path_list:
                                    if child['stId'] not in path_list:
                                        path_list.append(child['stId'])
                                        print(f"DBUG_GREAT: Greatgrandchild desired, child ID grabbed\t---\t"
                                              f"{path_list}")
        if child['type'] == 'Pathway' and child['stId'] in path_list:  # and child['stId'] in path_list
            print(f'DBUG3.1 --- Name: {child["name"]} --- stID: {child["stId"]} --- Type: Pathway --- Query: '
                  f'eventsHierarchy subtree')
            # rxn_dict[child['name']] = {}
            # get_path_data(child, rxn_dict[child['name']], path_list)
            get_path_data(child, path_list)
        elif (child['type'] == 'Reaction' or child['type'] == 'BlackBoxEvent') and rxn_flag is False:
            term_path_adapter(sub_dict)
            rxn_flag = True
    return

    for child in sub_dict['children']:
        if child['stId'] in path_list:



'''else:
    for child in sub_dict['hasEvent']:
        if child['schemaClass'] == 'Pathway':
            print(f'DBUG3.2: {child}')
            rxn_dict[child['displayName']] = {}
            path_url = f'/data/query/{child["stId"]}'
            full_url = url_base + path_url
            response = requests.get(full_url, headers=headers).json()
            print(f"DBUG4.2: QueryNeeded: {response}")
            get_path_data(response, rxn_dict[child['displayName']], 0)
        elif child['schemaClass'] == 'Reaction' or child['schemaClass'] == 'BlackBoxEvent':
            print(f'DBUG5.1: {child}')
            rxn_dict[child['displayName']] = {}
            get_parts_data(child, rxn_dict[child['displayName']])
    # else if child.type == "BlackBoxEvent":
    #    child['children'] = {}
    #    get_ortho_data(child)'''


# Base function to which all other functions are subordinate
def get_hier_data(entryId):
    path_url = '/data/eventsHierarchy/4530'
    full_url = url_base+path_url
    response = requests.get(full_url, headers=headers).json()
    print(f'DBUG1.1 --- Species: Oryza sativa --- taxID: 4530 --- Query: eventsHierarchy')
    # print(type(response))
    base_dict = {}
    for TopLevel in response:
        if TopLevel['name'] == "Metabolism and regulation":
            base_dict = TopLevel
            break
    print(f'DBUG2.1 --- Name: {base_dict["name"]} --- stID: {base_dict["stId"]} --- Query: eventsHierarchy subtree')
    rxn_dict = {}
    get_path_data(base_dict, entryId)
    with open(outDict, 'rb+') as filehandle:
        filehandle.seek(-1, os.SEEK_END)
        filehandle.truncate()
    outfile = open(outDict, "a")
    outfile.write("\n}")
    outfile.close()
    print(f'DBUG10.1 --- {rxn_dict}')
    return rxn_dict


# path_list = ['R-OSA-2744343', 'R-OSA-5655122', 'R-OSA-1119330', "R-OSA-1119319", "R-OSA-1119263",
#              "R-OSA-1119622", "R-OSA-1119289", "R-OSA-1119553", "R-OSA-1119354", "R-OSA-1119281",
#              "R-OSA-1119528", "R-OSA-1119567", "R-OSA-1119445", "R-OSA-1119297", "R-OSA-1119444"]
'''def get_hier_data(entryId):
    path_url = f'/data/query/{entryId}'
    full_url = url_base + path_url
    response = requests.get(full_url, headers=headers).json()
    print(f'DBUG1: {response}')
    rxn_dict = {}
    product_types = ['CatalystActivity', 'DefinedSet', 'Complex', 'EntityWithAccessionedEvent']
    if response['schemaClass'] == 'Pathway' or response['schemaClass'] == 'TopLevelPathway':
        print(f'DBUG2: {response["schemaClass"]}')
        get_path_data(response, rxn_dict, 0)
    elif response['schemaClass'] == 'Reaction' or response['schemaClass'] == 'BlackBoxEvent':
        get_parts_data(response, rxn_dict)
    elif response['schemaClass'] in product_types:
        get_product_data(response['dbId'], rxn_dict)
    return rxn_dict'''


def prettyPrint(rxn_dict, depth, outfile):
    if depth == 0:
        outfile.write(f"Pathway\tReaction\tUniProt_ID\tMSU?\tRAP?\tSpecies\tOrthologs\n")
    if 'species' in rxn_dict:
        extra = 0
        if 'MSU' in rxn_dict:
            outfile.write(f"{rxn_dict['MSU']}\t")
            extra += 4
        else:
            outfile.write(f"No MSU\t\t\t")
            extra += 4
        if 'RAP' in rxn_dict:
            outfile.write(f"{rxn_dict['RAP']}\t")
            extra += 4
        else:
            outfile.write(f"No RAP\t\t\t")
            extra +=4
        prettyPrint(rxn_dict['species'], depth + extra, outfile)
        # outfile.write(f'\n')
        # print(depth)
        if not len(rxn_dict['species']) > 0 and len(rxn_dict) > 0:
            outfile.write(f'\n')
    else:
        rxn_keys = [*rxn_dict]
        for element in rxn_dict:
            # print(f"{element} Depth: {depth}")
            if element != rxn_keys[0]:
                for tabs in range(0, depth):
                    outfile.write(f"\t")
            outfile.write(f"{element}\t")
            extra = math.floor(len(element)/4) + 1
            # print(f"{element} Extra Depth added: {extra}")
            if type(rxn_dict) == dict:
                prettyPrint(rxn_dict[element], depth + extra, outfile)
            else:
                outfile.write(f'\n')
            # else:
            #   outfile.write(f'\n')
                # print(bigstring)


if __name__ == '__main__':
    sys.argv = ['orthology_data_grabber.py', 'ortho_RPP_inter.json', 'ortho_DF_inter.csv', 'R-OSA-1119263']
    outDict = sys.argv[1]
    outFrame = sys.argv[2]
    pathways = sys.argv[3:]

    df_handle = open(outDict, "w")
    df_handle.write("{")
    df_handle.close()

    df_handle = open(outFrame, "w")
    df_handle.write("")
    df_handle.close()

    start_time = time.time()
    comp_struct = get_hier_data(pathways)
    print(f'finally: {comp_struct}')
    time1 = time.time() - start_time

    print("--- %s seconds ---" % time1)

    # input("Press Enter to continue...")
    # start_time1 = time.time()
    # df_handle = open(fileout, "w")
    # prettyPrint(comp_struct, 0, df_handle)
    # print(stringy)
    # time2 = time.time() - start_time1
    # print("--- %s seconds ---" % time2)
    # print("--- Time ratio: %f ---" % (time1/time2))

    exit(1)


'''def get_species_data():
    path_url = 'data/species/all'
    full_url = url_base + path_url
    response = requests.get(full_url, headers=headers).json
    return response[0]'''

# params = ('species','4530')

'''path_url = '/data/pathways/top/4530'
start_here = url_base+path_url
response = requests.get(start_here,
                        headers=headers,
                        #params=params,
                        )
for pathway in response:
    if pathway.displayName == "Metabolism and regulation":
        metabID = pathway.stID'''
