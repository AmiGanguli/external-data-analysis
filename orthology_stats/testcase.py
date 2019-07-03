#!
import requests
import re
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
    UniProtId = ewas['referenceEntity']['identifier']
    rxn_dict[UniProtId] = {}
    print(f'8: {ewas}')
    for name in ewas['referenceEntity']['geneName']:
        match = re.match("OS..G........", name)
        if match:
            rxn_dict[UniProtId]["RAP"] = name[0:12]
        match1 = re.match("LOC_OS..G.....", name)
        if match1:
            rxn_dict[UniProtId]['MSU'] = name
    rxn_dict[UniProtId]['species'] = {}
    if 'inferredTo' in ewas:
        for ortholog in ewas['inferredTo']:
            print(f'9: {ortholog}')
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
    print(f'7: {new_response}')
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
            print(f'6: {party}')
            get_product_data(party['peDbId'], rxn_dict)
    return


# Dives into pathways and constructs pathway hierarchy
def get_path_data(sub_dict, rxn_dict):
    for child in sub_dict['children']:
        print(f'3.5: {child}')
        if child['type'] == 'Pathway' and (child['stId'] == 'R-OSA-1119263' or child['stId'] == 'R-OSA-2744343' or child['stId'] == 'R-OSA-5655122'):
            print(f'4: {child}')
            rxn_dict[child['name']] = {}
            get_path_data(child, rxn_dict[child['name']])
        elif child['type'] == 'Reaction' or child['type'] == 'BlackBoxEvent':
            print(f'5: {child}')
            rxn_dict[child['name']] = {}
            get_parts_data(child, rxn_dict[child['name']])
        # else if child.type == "BlackBoxEvent":
        #    child['children'] = {}
        #    get_ortho_data(child)
    return


# Base function to which all other functions are subordinate
def get_hier_data():
    path_url = '/data/eventsHierarchy/4530'
    full_url = url_base+path_url
    response = requests.get(full_url, headers=headers).json()
    print(f'1: {response}')
    # print(type(response))

    base_dict = {}
    for TopLevel in response:
        if TopLevel['name'] == "Metabolism and regulation":
            base_dict = TopLevel
            break
    print(f'2: {base_dict}')
    rxn_dict = {}
    get_path_data(base_dict, rxn_dict)
    print(f'3: {rxn_dict}')
    return rxn_dict


comp_struct = get_hier_data()
print(f'finally: {comp_struct}')

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









