import asyncio
import pandas as pd
import time
from .schema import EventsHierarchy


class Data:
    def __init__(self, connection, saved_orthologs):
        self.connection = connection
        self.events_hierarchy = {}
        self.species_list = None

        # FIXME: This is now used for leaf pathway participants as well.
        # As the code is cleaned-up this should be made clearer.
        self.reaction_participants = {}
        self.reaction_orthologs = saved_orthologs
        self.reference_entities = None
        self.product_data = {}
        self.product_data_pending = set()

    async def eventsHierarchy(self, tax_id):
        if tax_id not in self.events_hierarchy:
            raw_events = await self.connection.getEventsHierarchy(tax_id)
            self.events_hierarchy[tax_id] = EventsHierarchy(raw_events)
        return self.events_hierarchy[tax_id]

    async def species(self):
        if self.species_list is None:
            species_raw = await self.connection.getSpecies()
            species_list = []
            for s in species_raw:
                species_list.append(
                    (s['dbId'], s['displayName'], s['name'], s['taxId'], s['abbreviation']))
            self.species_list = pd.DataFrame(
                data=species_list,
                columns=['id', 'display_name',
                         'names', 'tax_id', 'abbreviation']
            )
        return self.species_list

    async def participants(self, reaction_id):
        if reaction_id not in self.reaction_participants:
            self.reaction_participants[reaction_id] = await self.connection.getParticipantsPhysicalEntities(reaction_id)
            #self.reaction_participants[reaction_id] = await self.connection.getParticipants(reaction_id)
        print(f'Return participants for {reaction_id}')
        return self.reaction_participants[reaction_id]

    # Fetch a multiple of 20 ids.
    #

    async def product_data_fetch_block(self, priority_ids, fetch_incomplete_block=False):
        # We want an ordered list with your priority_ids at the start.
        self.product_data_pending -= priority_ids
        pending_ids = list(priority_ids) + list(self.product_data_pending)
        if len(pending_ids) == 0:
            return
        leave_for_later = len(pending_ids) % 20
        fetch_now = len(pending_ids) - leave_for_later
        if fetch_incomplete_block and fetch_now == 0:
            fetch_now = leave_for_later
            leave_for_later = 0

        ids_to_fetch = set(pending_ids[0:fetch_now])
        self.product_data_pending = set(pending_ids[fetch_now:])

        if fetch_now == 0:
            return

        ts = time.time()
        print(f'{ts} Batched fetch: {fetch_incomplete_block} {len(ids_to_fetch)}')
        async for record in self.connection.getProductDataMultiple(ids_to_fetch):
            self.product_data[record['dbId']] = record
        print(
            f'{ts} Batched fetch: {fetch_incomplete_block} {len(ids_to_fetch)} took {time.time() - ts}')

   
    async def productData(self, ids):
        missing_ids = ids.copy()
        for iteration in ['fetch bulk', 'fetch rest', 'done']:
            found_ids = set()
            for id in missing_ids:
                if id in self.product_data:
                    yield self.product_data[id]
                    found_ids.add(id)
            missing_ids -= found_ids
            if len(missing_ids) == 0:
                break
            while self.connection.willBlock():
                await asyncio.sleep(1)                
            if iteration == 'fetch bulk':
                await self.product_data_fetch_block(missing_ids)
            elif iteration == 'fetch rest':
                await self.product_data_fetch_block(missing_ids, fetch_incomplete_block=True)
            elif iteration == 'done':
                print(
                    'Error: somehow we reached the end of "done" while still missing ids {missing_ids}')
                exit(1)

    # Descend into reaction participants.
    # Given a participant's data record, check if it's a defined set.
    # if not, yield the data for the partipant.  If it is a defined
    # set, recurse into its members.
    #
    async def expandDefinedSets(self, participant_records):
        defined_sets = set()
        for participant in participant_records:
            # For some reason some compontent members are integers.
            # Need to investigate what's up with the API here.
            if type(participant) is int:
                continue
            if participant['schemaClass'] in ['DefinedSet', 'Complex']:
                defined_sets.add(participant['dbId'])
            else:
                yield participant
        if len(defined_sets) == 0:
            return
        # We've got a list of defined sets and complexes.  We need to recurse into
        # each one.
        async for record in self.productData(defined_sets):
            if 'hasMember' in record:
                members = self.expandDefinedSets(record['hasMember'])
            elif 'hasComponent' in record:
                members = self.expandDefinedSets(record['hasComponent'])
            else:
                print(
                    f'Strange, missing hasMember or hasComponent {participant["dbId"]}')
                continue
            async for member in members:
                yield member

    async def reactionOrthologs(self, parent, reaction):
        if reaction.stId in self.reaction_orthologs:
            return self.reaction_orthologs[reaction.stId]
        participants = await self.participants(reaction.stId)
        participants = self.expandDefinedSets(participants)
        results = []
        async for participant in participants:
            participant_id = participant['dbId']
            if participant['schemaClass'] == 'SimpleEntity':
                continue
            elif participant['schemaClass'] != 'EntityWithAccessionedSequence':
                print(
                    f"Don't know how to handle reaction participent of {participant['schemaClass']}")
                continue
            async for participant_data in self.productData(set([participant_id])):
                if 'referenceEntity' not in participant_data:
                    print(f'No referenceEntity for {participant_id}')
                    continue
                reference_entity = participant_data['referenceEntity']
                if 'databaseName' not in reference_entity:
                    print(f'No databasename for {participant_id}')
                    continue
                database_name = reference_entity['databaseName']
                if database_name != 'UniProt':
                    print(
                        f'Unexpected databaseName for {participant_id} is {database_name}')
                    continue
                uniprot_id = reference_entity['identifier']
                rap_id = reference_entity['geneName'][0]
                species_genes = {}
                if 'inferredTo' not in participant_data:
                    print(f'No orthologs in {participant_id}')
                else:
                    orthologs = self.expandDefinedSets(participant_data['inferredTo'])
                    async for ortholog in orthologs:
                        if ortholog['schemaClass'] == 'EntityWithAccessionedSequence':
                            species_name = ortholog['speciesName']
                            gene_name = ortholog['name'][0]
                            if species_name not in species_genes:
                                species_genes[species_name] = set()
                            species_genes[species_name].add(gene_name)
                results.append(
                    (parent, reaction, uniprot_id, rap_id, species_genes))
        self.reaction_orthologs[reaction.stId] = results
        return results

