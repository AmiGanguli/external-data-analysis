#!
import json
import requests
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

def get_product_data(sub_dict, rxn_dict):
    new_path = ('/data/participants/' + sub_dict.stId)
    new_full = url_base + new_path
    new_response = requests.get(new_full, headers=headers).json
    for party in new_response:
        if party.schemaClass != "SimpleEntity":
        #for entity in party.refEntities:
                #tggrgdrg
        if entity.schemaClass == "ReferenceGeneProduct":
            sub_dict[party.peDbId][entity.dbID] = {}
            sub_dict[party.peDbId][entity.dbID][UniProt] = entity.identifier
            for name in entity.geneName:
                if name.startswith("OS"):
                    sub_dict[party.peDbId][entity.dbID][RAP] = name
                if name.startswith("LOC"):
                    sub_dict[party.peDbId][entity.dbID][MSU] = name
            get_ortho_data()
    return

def get_ortho_data(event, species_dict):
    #for species in species_dict:
        id = species_dict.dbId
        new_path = ('/data/orthology/' + event.stId + '/species/' + id)
        new_full = url_base + new_path
        new_response = requests.get(new_full, headers=headers).json
        if new_response =
    return

def get_path_data(rxn_dict, species_dict):
    for child in rxn_dict.children:
        if child.type == "Pathway":
            get path_data(child, species_dict)
        else if child.type == "Reaction":
            child.children = {}
            get_ortho_data(child, species_dict)
        else if child.type == "BlackBoxEvent":
            child.children = {}
            get_ortho_data(child, species_dict)

def get_species_data():
    path_url = 'data/species/all'
    full_url = url_base + path_url
    response = requests.get(full_url, headers=headers).json
    return response[0]

#params = ('species','4530')

'''path_url = '/data/pathways/top/4530'
start_here = url_base+path_url
response = requests.get(start_here,
                        headers=headers,
                        #params=params,
                        )
for pathway in response:
    if pathway.displayName == "Metabolism and regulation":
        metabID = pathway.stID'''

#path_url = '/data/pathway/R-OSA-2744345/containedEvents'
path_url = '/data/eventsHierarchy/4530'
full_url = url_base+path_url
response = requests.get(full_url,headers=headers).json

rxn_dict = {}
for TopLevel in response:
    if TopLevel.stId == "R-OSA-2744345":
        rxn_dict = TopLevel
        break
#rxn_dict = {}
get_path_data(rxn_dict)

species_dict = get_species_data()







