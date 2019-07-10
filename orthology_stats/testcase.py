#!
import requests
import re
import sys
import math
import time
url_base = 'https://plantreactome.gramene.org/ContentService'
headers = {'accept': 'application/json'}

'''def get_path_data(rxn_dict):
    new_path = ('/data/pathway/' + rxn_dict.stID + '/containedEvents')
    new_full = url_base + new_path
    new_response = requests.get(new_full, headers=headers).json
    rxn_only = 1
    for sub_event in new_response:
        if sub_event.schemaClass == "Pathway":
            #get_path_data(sub_event.stID, rxn_dict)
            rxn_only = 0
    if rxn_only == 1:
        if sub_event.stID not in rxn_dict:
            rxn_dict[sub_event.stID] = {}
            for rxn in rxn_dict[sub_event.stID]:
                get_rxn_data(rxn.stID, rxn_dict[sub_event.stID])
    return'''
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

# path_url = '/data/pathway/R-OSA-2744345/containedEvents'


def get_multi_product_data(setId,rxn_list):
    new_path = (f'/data/query/{setId}')
    new_full = url_base + new_path
    new_response = requests.get(new_full, headers=headers).json()
    for member in new_response['hasMember']:
        rxn_list.append(member['name'][0])
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


def get_ortho_data(ewas, rxn_dict):
    if 'identifier' in ewas['referenceEntity']:
        UniProtId = ewas['referenceEntity']['identifier']
    else:
        UniProtId = ewas['referenceEntity']['secondaryIdentifier'][0]
    rxn_dict[UniProtId] = {}
    # print(f'DBUG8: {ewas}')
    if 'geneName' in ewas['referenceEntity']:
        for name in ewas['referenceEntity']['geneName']:
            match = re.match("OS..G........", name)
            if match:
                rxn_dict[UniProtId]["RAP"] = name[0:12]
            match1 = re.match("LOC_OS..G.....", name)
            if match1:
                rxn_dict[UniProtId]['MSU'] = name
    if 'inferredTo' in ewas:
        rxn_dict[UniProtId]['species'] = {}
        for ortholog in ewas['inferredTo']:
            # print(f'DBUG9: {ortholog}')
            rxn_dict[UniProtId]['species'][ortholog['speciesName']] = []
            if ortholog['schemaClass'] == "EntityWithAccessionedSequence":
                rxn_dict[UniProtId]['species'][ortholog['speciesName']].append(ortholog['name'][0])
            if ortholog['schemaClass'] == 'DefinedSet':
                get_multi_product_data(ortholog['stId'], rxn_dict[UniProtId]['species'][ortholog['speciesName']])
    return


def get_product_data(entityId, rxn_dict):
    new_path = (f'/data/query/{entityId}')
    new_full = url_base + new_path
    new_response = requests.get(new_full, headers=headers).json()
    print(f'DBUG7: {new_response}')
    if new_response['schemaClass'] == "CatalystActivity":
        get_product_data(new_response['physicalEntity']['stId'], rxn_dict)
        print('-> Catalyst')
    elif new_response['schemaClass'] == "DefinedSet":
        print('-> Set')
        for member in new_response['hasMember']:
            print(f'--->{member["schemaClass"]}')
            if member['schemaClass'] != "SimpleEntity":
                get_product_data(member['stId'], rxn_dict)
    elif new_response['schemaClass'] == 'Complex':
        print('->Complex')
        for component in new_response['hasComponent']:
            if component['schemaClass'] != "SimpleEntity":
                print(f'--->{component["schemaClass"]}')
                get_product_data(component['stId'], rxn_dict)
    elif new_response['schemaClass'] == "EntityWithAccessionedSequence":
        print('->EWAS')
        get_ortho_data(new_response, rxn_dict)
    return


# Something of a wrapper function to reach CatalystActivity
def get_parts_data(event, rxn_dict):
    new_path = ('/data/participants/' + event['stId'])
    new_full = url_base + new_path
    new_response = requests.get(new_full, headers=headers).json()
    for party in new_response:
        if party['schemaClass'] == "CatalystActivity":
            print(f'DBUG6: {party}')
            get_product_data(party['peDbId'], rxn_dict)
    return


# (child['stId'] == 'R-OSA-1119263' or
# Dives into pathways and constructs pathway hierarchy
def get_path_data(sub_dict, rxn_dict, path_list):
    if path_list:
        for child in sub_dict['children']:
            # print(f'3.5: {child}')
            if child['stId'] in path_list:
                if 'children' in child:
                    for grandkid in child['children']:
                        if grandkid['type'] == 'Pathway':
                            path_list.append(grandkid['stId'])
                            print(f"DBUG: Child reached, grandchildren IDs grabbed")
                            print(path_list)
            if child['stId'] not in path_list:
                if 'children' in child:
                    for grandkid in child['children']:
                        if grandkid['stId'] in path_list:
                            if child['stId'] not in path_list:
                                path_list.append(child['stId'])
                                print(f"DBUG: Grandchild desired, child ID grabbed")
                                print(path_list)
            if child['type'] == 'Pathway' and child['stId'] in path_list:  # and child['stId'] in path_list
                print(f'DBUG4: {child}')
                rxn_dict[child['name']] = {}
                get_path_data(child, rxn_dict[child['name']], path_list)
            elif child['type'] == 'Reaction' or child['type'] == 'BlackBoxEvent':
                print(f'DBUG5: {child}')
                rxn_dict[child['name']] = {}
                get_parts_data(child, rxn_dict[child['name']])
    else:
        for child in sub_dict['hasEvent']:
            if child['schemaClass'] == 'Pathway':
                print(f'DBUG4: {child}')
                rxn_dict[child['displayName']] = {}
                path_url = f'/data/query/{child["stId"]}'
                full_url = url_base + path_url
                response = requests.get(full_url, headers=headers).json()
                print(f"DBUG3.2: QueryNeeded: {response}")
                get_path_data(response, rxn_dict[child['displayName']], 0)
            elif child['schemaClass'] == 'Reaction' or child['schemaClass'] == 'BlackBoxEvent':
                print(f'DBUG5: {child}')
                rxn_dict[child['displayName']] = {}
                get_parts_data(child, rxn_dict[child['displayName']])
        # else if child.type == "BlackBoxEvent":
        #    child['children'] = {}
        #    get_ortho_data(child)
    return


# Base function to which all other functions are subordinate
def get_hier_data(entryId):
    path_url = '/data/eventsHierarchy/4530'
    full_url = url_base+path_url
    response = requests.get(full_url, headers=headers).json()
    print(f'DBUG1.1: {response}')
    # print(type(response))
    base_dict = {}
    for TopLevel in response:
        if TopLevel['name'] == "Metabolism and regulation":
            base_dict = TopLevel
            break
    print(f'DBUG2.1: {base_dict}')
    rxn_dict = {}
    get_path_data(base_dict, rxn_dict, entryId)
    print(f'DBUG3.1: {rxn_dict}')
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
        if len(rxn_dict['species']) > 0 and len(rxn_dict) > 0:
            # for tabs in range(0, depth):
            #    outfile.write(f"\t")
            print('to test this')
        else:
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
        rxn_keys = [*rxn_dict]
        for element in rxn_dict:
            print(f"{element} Depth: {depth}")
            outfile.write(f"{element}\t")
            extra = math.floor(len(element)/4) + 1
            print(f"{element} Extra Depth added: {extra}")
            if type(rxn_dict) == dict:
                usefulPrint(rxn_dict[element], path_depth-1, depth + extra, outfile)
            else:
                outfile.write(f'\n')


sys.argv = ['testcase.py', 'R-OSA-5655122']
fileout = "testcaseout.txt"
start_time = time.time()
comp_struct = get_hier_data([sys.argv[1]])
print(f'finally: {comp_struct}')
print("--- %s seconds ---" % (time.time() - start_time))
input("Press Enter to continue...")
start_time = time.time()
df_handle = open(fileout, "w")
prettyPrint(comp_struct, 0, df_handle)
# print(stringy)
print("--- %s seconds ---" % (time.time() - start_time))
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