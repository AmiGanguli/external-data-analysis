import asyncio
import pandas as pd
from .schema import EventsHierarchy


class Data:
    def __init__(self, connection):
        self.connection = connection
        self.events_hierarchy = {}
        self.species_list = None
        self.reaction_participants = {}
        self.reference_entities = None

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
        return self.reaction_participants[reaction_id]
