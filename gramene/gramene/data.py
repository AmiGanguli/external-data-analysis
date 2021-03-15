import asyncio
import pandas as pd
import time
from .schema import EventsHierarchy


class Data:
    def __init__(self, connection):
        self.connection = connection
        self.events_hierarchy = {}
        self.species_list = None
        self.reaction_participants = {}
        self.reference_entities = None
        self.product_data = {}
        self.product_data_in_flight = set()
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
        print(f'Return participants for {reaction_id}')
        return self.reaction_participants[reaction_id]

    async def product_data_fetch_next_20(self):
        if len(self.product_data_pending) == 0:
            return
        pending = list(self.product_data_pending)
        ids_to_fetch = pending[0:20]
        ts = time.time()
        print(f'{ts} Batched fetch: {len(ids_to_fetch)}')
        ids_to_fetch_set = set(ids_to_fetch)
        self.product_data_pending -= ids_to_fetch_set
        self.product_data_in_flight |= ids_to_fetch_set
        records = await self.connection.getProductDataMultiple(ids_to_fetch_set)
        print(f'{ts} Batched fetch: {len(ids_to_fetch)} took {time.time() - ts}')
        self.product_data_in_flight -= ids_to_fetch_set
        for record in records:
            self.product_data[record['dbId']] = record

    # Attempts to batch 20 product data calls together.
    #
    # The intent is for many of these to be running (sort of) simultaneously.  While
    # product_data_fetch_next_20 is blocked, other threads can be adding to
    # product_data_pending. This increases the likelihood of sending "full" requests
    # without adding any logic to coordinate the threads beyond the shared pending list.
    #
    # FIXME: We could get more parrallelism out of this by getting multiple blocks of
    # 20 at the same time.  Not sure if it's worth it, but something to try.
    #
    async def productData(self, ids):
        missing_ids = ids.copy()
        done = False
        while not done:
            found_ids = set()
            required_ids = set()
            requests_in_flight = False
            for id in missing_ids:
                if id in self.product_data:
                    yield self.product_data[id]
                    found_ids.add(id)
                elif id in self.product_data_in_flight:
                    requests_in_flight = True
                elif id not in self.product_data_pending:
                    required_ids.add(id)
            missing_ids -= found_ids
            self.product_data_pending |= required_ids
            if len(missing_ids) > 0:
                if requests_in_flight:
                    await asyncio.sleep(0.5)
                elif len(self.product_data_pending) > 0:
                    print(
                        f'fetching next 20 with {len(missing_ids)} missing ids and {len(self.product_data_pending)} pending requests')
                    await self.product_data_fetch_next_20()
            else:
                done = True
