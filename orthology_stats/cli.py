import click
import os
import re
import asyncio
from gramene.client import Connection
from gramene.data import Data
import pandas as pd
from species_list import species_list
import aiohttp
import logging

version = '0.01'

# Default options
api_endpoint = 'https://plantreactomedev.gramene.org/ContentService'
output_directory = './'


class Reactome:
    def __init__(self, tax_id, top_level, api_endpoint, output_directory, use_saved_orthologs, log_file, log_level, file_prefix, show):
        self.tax_id = tax_id
        self.top_level = top_level
        self.api_endpoint = api_endpoint
        self.output_directory = output_directory
        self.use_saved_orthologs = use_saved_orthologs
        self.file_prefix = file_prefix
        self.show = show
        pd.set_option('display.max_rows', None)
        level = {
            'DEBUG': logging.DEBUG,
            'INFO': logging.INFO,
            'WARNING': logging.WARNING,
            'ERROR': logging.ERROR,
            'CRITICAL': logging.CRITICAL,
        }[log_level]
        if log_file is not None:
            logging.basicConfig(format='%(asctime)s %(message)s',
                                filename=log_file, level=level)
        else:
            logging.basicConfig(format='%(asctime)s %(message)s',
                                filename=log_file, level=level)


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
@click.option(
    '--use-saved-orthologs',
    envvar='REACTOME_SAVED_ORTHOLOGS_FILE',
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    help='Read orthologs from this file rather than fetching remotely.'
)
@click.option(
    '--log-file',
    envvar='REACTOME_LOG_FILE',
    type=click.Path(file_okay=True, writable=True),
    help='Log to this file.'
)
@click.option(
    '--log-level',
    envvar='REACTOME_LOG_LEVEL',
    type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR',
                       'CRITICAL'], case_sensitive=False),
    default='WARNING',
    help='Log level.'
)
@click.option('--file-prefix', default='', help='Prepend this string to output files.')
@click.option('--show/--no-show', default=True, help='Display the output or a summary of the output.')
@click.pass_context
def cli(ctx, tax_id, top_level, api_endpoint, output_directory, use_saved_orthologs, log_file, log_level, file_prefix, show):
    ctx.obj = Reactome(tax_id, top_level, api_endpoint,
                       output_directory, use_saved_orthologs, log_file, log_level, file_prefix, show)




async def species_(ctx):
    species_df = None
    async with aiohttp.ClientSession(headers=headers) as session:
        connection = Connection(session, ctx.obj.api_endpoint)
        data = Data(connection)
        species_df = await data.species()
    return species_df


@cli.command()
@click.pass_context
def species(ctx):
    species_df = asyncio.run(species_(ctx))
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

def run_with_connection(func):
    async def inner(*args, **kwargs):
        result = None
        ctx = args[0]
        connector = aiohttp.TCPConnector(limit_per_host=15)
        async with aiohttp.ClientSession(connector=connector, headers=headers) as session:
            connection = Connection(session, ctx.obj.api_endpoint, 15)
            data = Data(connection)
            kwargs['connection'] = connection
            kwargs['data'] = data
            result = await func(*args, **kwargs)
        return result
    def runner(*args, **kwargs):
        return asyncio.run(inner(*args, **kwargs))
    return runner

async def get_events(obj, data):
    events = await data.eventsHierarchy(obj.tax_id)
    if obj.top_level is not None:
        events = events[obj.top_level]
    return events

headers = {'accept': 'application/json', 'content-type': 'text/plain'}

@run_with_connection
async def eventtree_(ctx, connection, data):
    return await get_events(ctx.obj, data)


@cli.command()
@click.pass_context
def eventtree(ctx):
    events = eventtree_(ctx)
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


@run_with_connection
async def allOrthologs_(ctx, connection, data):
    events = None
    result = []
    events = await get_events(ctx.obj, data)

    orthologs = []
    for parent, reaction in events.all_reactions():
        orthologs.append(data.reactionOrthologs(parent, reaction))
    orthologs_results = await asyncio.gather(*orthologs)
    for ortholog_result in orthologs_results:
        for parent, reaction, uniprot_id, rap_id, species_genes in ortholog_result:
            result.append((parent.stId, reaction.stId, uniprot_id, rap_id, species_genes))
    return result


@cli.command()
@click.pass_context
def all_orthologs(ctx):
    result = allOrthologs_(ctx)

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
