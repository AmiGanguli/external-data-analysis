import click
import os
import re
import asyncio
from gramene.client import Connection
from gramene.data import Data
import pandas as pd
from species_list import species_list
import aiohttp

version = '0.01'

# Default options
api_endpoint = 'https://plantreactomedev.gramene.org/ContentService'
output_directory = './'


class Reactome:
    def __init__(self, tax_id, top_level, api_endpoint, output_directory, file_prefix, show, debug):
        self.tax_id = tax_id
        self.top_level = top_level
        self.api_endpoint = api_endpoint
        self.output_directory = output_directory
        self.file_prefix = file_prefix
        self.show = show
        self.debug = debug
        self.debug_message('Starting...')
        pd.set_option('display.max_rows', None)

    def debug_message(self, message):
        if self.debug:
            print(message)


@click.group(chain=True)
@click.option(
    '--tax-id',
    default=4530,
    type=click.INT,
    help='NCBI Tax ID for the reference organism.'
)
@click.option(
    '--top-level',
    type=click.STRING,
    help='Top level pathway name.'
)
@click.option(
    '--api-endpoint',
    envvar='REACTOME_API',
    default=api_endpoint,
    help='The url for the Reactome API.'
)
@click.option(
    '--output-directory',
    envvar='REACTOME_OUTPUT_DIRECTORY',
    type=click.Path(exists=True, file_okay=False, writable=True),
    help='The destintation for output files.'
)
@click.option('--file-prefix', default='', help='Prepend this string to output files.')
@click.option('--show/--no-show', default=True, help='Display the output or a summary of the output.')
@click.option('--debug/--no-debug', envvar='REACTOME_DEBUG', default=False)
@click.pass_context
def cli(ctx, tax_id, top_level, api_endpoint, output_directory, file_prefix, show, debug):
    ctx.obj = Reactome(tax_id, top_level, api_endpoint,
                       output_directory, file_prefix, show, debug)


async def get_events(obj, data):
    events = await data.eventsHierarchy(obj.tax_id)
    if obj.top_level is not None:
        events = events[obj.top_level]
    return events

headers = {'accept': 'application/json', 'content-type': 'text/plain'}


async def eventtree_(ctx):
    events = None
    async with aiohttp.ClientSession(headers=headers) as session:
        connection = Connection(session, ctx.obj.api_endpoint)
        data = Data(connection)
        events = await get_events(ctx.obj, data)
    return events


@cli.command()
@click.pass_context
def eventtree(ctx):
    events = asyncio.run(eventtree_(ctx))
    if events is None:
        print(f'No events at {ctx.obj.top_level}.')
        return
    if ctx.obj.show:
        for event, depth in events.walk():
            print(f'{depth * "  "} {event}')
    if ctx.obj.output_directory:
        out = events.to_data_frame()
        out.to_csv(
            path_or_buf=os.path.join(
                ctx.obj.output_directory, ctx.obj.file_prefix + 'event_hierarchy.csv'),
            mode='w',
            index_label='Row'
        )


@cli.command()
@click.pass_context
def species(ctx):
    species_df = asyncio.run(ctx.obj.data.species())
    if species_df is None:
        print(f'No species found.')
        return
    if ctx.obj.show:
        print(species_df)
    if ctx.obj.output_directory:
        species_df.to_csv(
            path_or_buf=os.path.join(
                ctx.obj.output_directory, ctx.obj.file_prefix + 'species.csv'),
            mode='w',
            index_label='Row'
        )


@cli.command()
@click.pass_context
def reactions(ctx):
    events = get_events(ctx.obj)
    if events is None:
        print(f'No events at {ctx.obj.top_level}.')
        return
    if ctx.obj.show:
        print(f'{"Parent pathway":20} {"Reaction":20}')
        print('-' * 20 + ' ' + '-'*20)
    for parent, reaction in events.all_reactions():
        if ctx.obj.show:
            print(f'{parent.stId:20} {reaction.stId:20}')


@cli.command()
@click.argument('id')
@click.pass_context
def reaction_participants(ctx, id):
    participants = asyncio.run(ctx.obj.data.participants(id))
    if ctx.obj.show:
        print(f'Reaction participants for {id}')
        print('-------------------------------')
        for participant in participants:
            print(participant)


def get_proteins(ctx):
    events = get_events(ctx.obj)
    if events is None:
        return
    for parent, reaction in events.all_reactions():
        participant_data = asyncio.run(
            ctx.obj.connection.getParticipantsReferenceEntities(reaction.stId))
        for p in participant_data:
            if 'databaseName' in p and p['databaseName'] == 'UniProt':
                yield parent.stId, reaction.stId, p['identifier'], str(p['dbId'])


def get_proteins_physical_entities(ctx):
    events = get_events(ctx.obj)
    if events is None:
        return
    for parent, reaction in events.all_reactions():
        participant_data = asyncio.run(
            ctx.obj.connection.getParticipantsPhysicalEntities(reaction.stId))
        for p in participant_data:
            if p['schemaClass'] == 'CatalystActivity':
                yield parent, reaction, p['schemaClass'], p['peDbId']
            elif p['schemaClass'] in ['DefinedSet', 'Complex', 'EntityWithAccessionedSequence']:
                yield parent, reaction, p['schemaClass'], p['stId']


@cli.command()
@click.pass_context
def malformed_protein_identifiers(ctx):
    if ctx.obj.show:
        print(f'{"Parent pathway":20} {"Reaction":20} {"UniProt ID":20}')
        print('-' * 20 + ' ' + '-'*20 + ' ' + '-'*20)
    for pathway, reaction, protein, protein_id in get_proteins(ctx):
        if len(protein) <= 6:
            continue
        if ctx.obj.show:
            print(f'{pathway:20} {reaction:20} {protein:20}')


@cli.command()
@click.pass_context
def list_proteins(ctx):
    if ctx.obj.show:
        print(f'{"Parent pathway":20} {"Reaction":20} {"UniProt ID":20}')
        print('-' * 20 + ' ' + '-'*20 + ' ' + '-'*20)
    data = []
    for pathway, reaction, protein, protein_id in get_proteins(ctx):
        data.append((pathway, reaction, protein))
        if ctx.obj.show:
            print(f'{pathway:20} {reaction:20} {protein:20} {protein_id:20}')

    if ctx.obj.output_directory:
        df = pd.DataFrame(data, columns=['Pathway', 'Reaction', 'Protein'])
        df.sort_values(by=['Pathway', 'Reaction', 'Protein'], inplace=True)
        df.to_csv(
            path_or_buf=os.path.join(
                ctx.obj.output_directory, ctx.obj.file_prefix + 'proteins.csv'),
            mode='w',
            index_label='Row'
        )


@cli.command()
@click.pass_context
def list_proteins_with_species(ctx):
    proteins = set()
    rows = []
    for pathway, reaction, protein, protein_id in get_proteins(ctx):
        proteins.add(protein_id)
        rows.append((pathway, reaction, protein))
    products = asyncio.run(
        ctx.obj.connection.getProductDataMultiple(proteins))

    # The key is species['dbId']-protein['identifier']
    # Value is a set of gene ids
    # if this works, we still need to orthologs
    genes = {}
    product_ids = set()
    species_ids = set()
    for product in products:
        if 'identifier' not in product:
            print('Missing identifier for ', product)
            continue
        if 'species' not in product:
            print('Missing species for ', product)
            continue
        if 'dbId' not in product:
            print('Missing dbId for species in ', product)
            continue
        if 'geneName' not in product:
            print('Missing geneName for ', product)
            continue
        identifier = product['identifier']
        species_id = str(product['species']['dbId'])
        product_ids |= {identifier, }
        species_ids |= {species_id, }
        gene = None
        for name in product['geneName']:
            match = re.match("OS..G.......*", name.upper())
            if match:
                gene = name
                break
        if gene is None:
            print(f'no gene found for {identifier}')
            continue
        key = species_id + '-' + identifier
        if key not in genes:
            genes[key] = ()
        genes[key] += (gene,)
    if ctx.obj.show:
        print(f'{"Parent pathway":20} {"Reaction":20} {"UniProt ID":20}')
        print('-' * 20 + ' ' + '-'*20 + ' ' + '-'*20)
        for row in rows:
            print(f'{row[0]:20} {row[1]:20} {row[2]:20}')
            for species_id in species_ids:
                key = species_id + '-' + row[2]
                if key in genes:
                    print('       ', key, ': ', genes[key])


# Descend into reaction participants.
# Given a participant's data record, check if it's a defined set.
# if not, yield the data for the partipant.  If it is a defined
# set, recurse into its members.
#


async def expand_defined_sets(ctx, connection, participant_records):
    for participant in participant_records:
        schema = participant['schemaClass']
        if schema == 'DefinedSet':
            set_data = await connection.getProductData(participant['dbId'])
            if 'hasMember' in set_data:
                members = expand_defined_sets(
                    ctx, connection, set_data['hasMember'])
                async for member in members:
                    yield member
            else:
                print(
                    f'Strange, defined set mising hasMember {participant["dbId"]}')
        else:
            yield participant


@cli.command()
@click.argument('id')
@click.pass_context
def reaction_participants_orthologs(ctx, id):
    participants = asyncio.run(ctx.obj.data.participants(id))
    result = []
    for participant in expand_defined_sets(ctx, participants):
        participant_id = participant['dbId']
        if participant['schemaClass'] == 'SimpleEntity':
            continue
        elif participant['schemaClass'] != 'EntityWithAccessionedSequence':
            print(
                f"Don't know how to handle reaction participent of {participant['schemaClass']}")
            continue
        participant_data = asyncio.run(
            ctx.obj.connection.getProductData(participant_id))
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
            for ortholog in expand_defined_sets(ctx, participant_data['inferredTo']):
                if ortholog['schemaClass'] == 'EntityWithAccessionedSequence':
                    species_name = ortholog['speciesName']
                    gene_name = ortholog['name'][0]
                    if species_name not in species_genes:
                        species_genes[species_name] = set()
                    species_genes[species_name].add(gene_name)
        result.append((uniprot_id, rap_id, species_genes))
    print(f'{"UniProt ID":20} {"RAP ID":20} {"Species":20}')
    for row in result:
        print(f'{row[0]:20} {row[1]:20}')
        for species_name, genes in row[2].items():
            print(' ' * 42 + species_name)
            for gene in genes:
                print(' ' * 50 + gene)

# Fixme: refactor this to reuse code from above.  Right now it's cut-and-paste.
#


def reactions(ctx):
    events = get_events(ctx.obj)
    if events is None:
        print(f'No events at {ctx.obj.top_level}.')
        return
    if ctx.obj.show:
        print(f'{"Parent pathway":20} {"Reaction":20}')
        print('-' * 20 + ' ' + '-'*20)
    for parent, reaction in events.all_reactions():
        if ctx.obj.show:
            print(f'{parent.stId:20} {reaction.stId:20}')


async def all_orthologs_(ctx):
    events = None
    result = []
    async with aiohttp.ClientSession(headers=headers) as session:
        connection = Connection(session, ctx.obj.api_endpoint)
        data = Data(connection)
        events = await get_events(ctx.obj, data)

        for parent, reaction in events.all_reactions():
            participants = await data.participants(reaction.stId)
            participants = expand_defined_sets(ctx, connection, participants)
            async for participant in participants:
                participant_id = participant['dbId']
                if participant['schemaClass'] == 'SimpleEntity':
                    continue
                elif participant['schemaClass'] != 'EntityWithAccessionedSequence':
                    print(
                        f"Don't know how to handle reaction participent of {participant['schemaClass']}")
                    continue
                participant_data = await connection.getProductData(participant_id)
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
                    orthologs = expand_defined_sets(
                        ctx, connection, participant_data['inferredTo'])
                    async for ortholog in orthologs:
                        if ortholog['schemaClass'] == 'EntityWithAccessionedSequence':
                            species_name = ortholog['speciesName']
                            gene_name = ortholog['name'][0]
                            if species_name not in species_genes:
                                species_genes[species_name] = set()
                            species_genes[species_name].add(gene_name)
                result.append(
                    (parent.stId, reaction.stId, uniprot_id, rap_id, species_genes))
    return result


@cli.command()
@click.pass_context
def all_orthologs(ctx):
    result = asyncio.run(all_orthologs_(ctx))

    print(f'{"Pathway ID":20} {"Reaction ID":20} {"UniProt ID":20} {"RAP ID":20} {"Species":20}')

    df = pd.DataFrame(columns=['Pathway', 'Reaction',
                               'Protein', 'RAP ID'] + species_list)

    for row in result:
        print(f'{row[0]:20} {row[1]:20} {row[2]:20} {row[3]:20}')
        new_row = {
            'Pathway': row[0],
            'Reaction': row[1],
            'Protein': row[2],
            'RAP ID': row[3],
        }
        for species_name, genes in row[4].items():
            if species_name in species_list:
                new_row[species_name] = '|'.join(genes)
            print(' ' * 42 + species_name)
            for gene in genes:
                print(' ' * 50 + gene)
        # FIXME: inefficient
        df = df.append(new_row, ignore_index=True)

    if ctx.obj.output_directory:
        df.sort_values(by=['Pathway', 'Reaction', 'Protein'], inplace=True)
        df.to_csv(
            path_or_buf=os.path.join(
                ctx.obj.output_directory, ctx.obj.file_prefix + 'orthologs.csv'),
            mode='w'
        )


if __name__ == '__main__':
    cli()
