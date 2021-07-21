import click
import os
import re
import asyncio
from gramene.client import Connection
from gramene.data import Data
import pandas as pd
import numpy as np
from species_list import species_list
import aiohttp
import logging

version = '0.01'

# Default options
api_endpoint = 'https://plantreactomedev.gramene.org/ContentService'
output_directory = './'


class Reactome:
    def __init__(self, tax_id, top_level, api_endpoint, output_directory, saved_orthologs, log_file, log_level, file_prefix, show):
        self.tax_id = tax_id
        self.top_level = top_level
        self.api_endpoint = api_endpoint
        self.output_directory = output_directory
        self.saved_orthologs = saved_orthologs
        self.file_prefix = file_prefix
        self.show = show

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

    # Global options
    pd.set_option('display.max_rows', None)

    saved_orthologs = {}
    if use_saved_orthologs is not None:
        saved_orthologs_df = pd.read_csv(use_saved_orthologs)
        species_names = list(saved_orthologs_df.columns)
        done = False
        while not done:
            if species_names[0] == 'RAP ID':
                done = True
            species_names = species_names[1:]
        for row in saved_orthologs_df.to_dict(orient='records'):
            parent = row['Pathway']
            reaction = row['Reaction']
            uniprot_id = row['Protein']
            rap_id = row['RAP ID']
            species_genes = {}
            for species_name in species_names:
                if pd.isna(row[species_name]):
                    continue
                species_genes[species_name] = row[species_name].split('|')
            if reaction not in saved_orthologs_df:
                saved_orthologs[reaction] = []
            saved_orthologs[reaction].append((parent, reaction, uniprot_id, rap_id, species_genes))

    # Create context object that will be passed to sub-commands
    ctx.obj = Reactome(tax_id, top_level, api_endpoint,
                       output_directory, saved_orthologs, log_file, log_level, file_prefix, show)


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
    events = eventtree_(ctx)
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
            data = Data(connection, ctx.obj.saved_orthologs)
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
        for event, depth, tree, child_prefix in events.walk():
            print(f'{tree}{event}')
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
    #result = []
    events = await get_events(ctx.obj, data)

    orthologs = []
    for parent, depth, tree, children_prefix in events.walk():
        for reaction in parent.reactions:
            orthologs.append(data.reactionOrthologs(parent, reaction))
        for reaction in parent.black_box_events:
            orthologs.append(data.reactionOrthologs(parent, reaction))
    orthologs_results = await asyncio.gather(*orthologs)

    ortholog_dict = {}
    for ortholog_result in orthologs_results:
        for parent, reaction, uniprot_id, protein_dbid, rap_id, species_genes in ortholog_result:
            if parent.stId not in ortholog_dict:
                ortholog_dict[parent.stId] = {}
            if reaction.stId not in ortholog_dict[parent.stId]:
                ortholog_dict[parent.stId][reaction.stId] = {}
            ortholog_dict[parent.stId][reaction.stId][uniprot_id] = (
                protein_dbid,
                rap_id,
                species_genes
            )

    return events, ortholog_dict


@cli.command()
@click.option('--count/--no-count', default=False)
@click.option('--show-empty/--no-show-empty', default=True)
@click.option('--tree/--no-tree', default=True)
@click.option('--proteins/--no-proteins', default=True)
@click.pass_context
def all_orthologs(ctx, count, show_empty, tree, proteins):
    events, ortholog_dict = allOrthologs_(ctx)

    print('Saving result')
    max_depth = events.max_depth()
    pathway_titles = [f'Pathway {i}' for i in range(max_depth+1)]

    df = pd.DataFrame(columns=['Tree'] + pathway_titles + ['Reaction',
            'Protein', 'Protein dbId', 'RAP ID', 'Orthologs Sum', 'Species Count'] + species_list)

    pathway_counts = {}
    reaction_counts = {}
    for parent, depth, tree, children_prefix, full_path in events.walk(full_path=[]):
        pathway_row = {f'Pathway {i}':full_path[i].stId for i in range(len(full_path))}
        pathway_row['Tree'] = tree + parent.name
        df = df.append(pathway_row, ignore_index=True)
        if parent.stId not in ortholog_dict:
            continue
        for reaction_id, products in ortholog_dict[parent.stId].items():
            reaction_row = pathway_row.copy()
            reaction_row['Tree'] = children_prefix + '   ' + reaction_id
            reaction_row['Reaction'] = reaction_id
            df = df.append(reaction_row, ignore_index=True);
            if not reaction_id in reaction_counts:
                reaction_counts[reaction_id] = {'species':0,'orthologs':0}
            for uniprot_id, orthologs in products.items():
                protein_row = reaction_row.copy()
                protein_row.update({
                    'Tree': children_prefix + '       ' + uniprot_id,
                    'Protein': uniprot_id,
                    'Protein dbId': orthologs[0],
                    'RAP ID': orthologs[1],
                    'Orthologs Sum': 0,
                    'Species Count': 0,
                })
                for species_name, genes in orthologs[2].items():
                    if species_name not in species_list:
                        continue
                    gene_count = len(genes)
                    protein_row['Species Count'] += 1
                    protein_row['Orthologs Sum'] += gene_count
                    if count:
                        protein_row[species_name] = gene_count
                    else:
                        protein_row[species_name] = '|'.join(genes)
                for path_segment in full_path:
                    segment_id = path_segment.stId
                    if not segment_id in pathway_counts:
                        pathway_counts[segment_id] = {'species':0,'orthologs':0}
                    pathway_counts[segment_id]['species'] += protein_row['Species Count']
                    pathway_counts[segment_id]['orthologs'] += protein_row['Orthologs Sum']
                    print('segment_id', segment_id)
                reaction_counts[reaction_id]['species'] += protein_row['Species Count']
                reaction_counts[reaction_id]['orthologs'] += protein_row['Orthologs Sum']
                if not show_empty and protein_row['Species Count'] == 0:
                    continue
                df = df.append(protein_row, ignore_index=True)

    df.drop_duplicates(keep=False, inplace=True)
    for index_label, row in df.iterrows():
        #pathway_titles = [f'Pathway {i}' for i in range(max_depth+1)]
        #pathway_row = {f'Pathway {i}':full_path[i].stId for i in range(len(full_path))}
        #print(type(row['Protein']).__name__, row['Protein'])
        if type(row['Reaction']) == str:
            if type(row['Protein']) == str:
                continue
            df.at[index_label, 'Species Count'] = reaction_counts[row['Reaction']]['species']
            df.at[index_label, 'Orthologs Sum'] = reaction_counts[row['Reaction']]['orthologs']
            continue
        #print('put in value for row', row)
        for i in range(max_depth, -1, -1):
            parent = row[f'Pathway {i}']
            if type(parent) == str:
                print('parent', parent)
                if parent not in pathway_counts:
                    print('missing parent from pathway counts', parent)
                else:
                    df.at[index_label, 'Species Count'] = pathway_counts[parent]['species']
                    df.at[index_label, 'Orthologs Sum'] = pathway_counts[parent]['orthologs']
                break

    if ctx.obj.output_directory:
        if count:
            file_suffix = 'orthologs-count'
        else:
            file_suffix = 'orthologs'
        if not proteins:
            df = df[df['Protein'].isna()]
        df.to_csv(
            path_or_buf=os.path.join(
                ctx.obj.output_directory, ctx.obj.file_prefix + file_suffix + '.csv'),
            index=False,
            mode='w'
        )

@run_with_connection
async def pathway_participants_(ctx, event_id, connection, data):
    result = {}
    participants = await data.participants(event_id)
    participants = data.expandDefinedSets(participants)
    async for participant in participants:
        #print(participant)
        if participant['className'] not in result:
            result[participant['className']] = set()
        result[participant['className']].add(participant['stId'])
    return result
        
@cli.command()
@click.argument('id')
@click.pass_context
def pathway_participants(ctx, id):
    classes = pathway_participants_(ctx, id)
    for class_name, items in classes.items():
        print(f'{class_name} {len(items)}')

    print(classes)

@cli.command()
@click.pass_context
def pathway_nodes_with_reactions(ctx):
    events = eventtree_(ctx)
    if events is None:
        print(f'No events at {ctx.obj.top_level}.')
        return

    columns = {'Pathway': [], 'tree': [], 'Reaction count': [], 'Protein count': [], 'DNA Sequence count': [], 'RNA Sequence count': []}
    trees = []
    names = []
    max_label = 0
    for event, depth, tree in events.walk():
        max_label = max(max_label, len(tree) + len(event.name))
        stId = '/'
        if event.name != 'TopLevel':
            stId = event.stId
            participants = pathway_participants_(ctx, stId)
            if 'Protein' in participants:
                columns['Protein count'].append(len(participants['Protein']))
            else:
                columns['Protein count'].append(0)
            if 'DNA Sequence' in participants:
                columns['DNA Sequence count'].append(len(participants['DNA Sequence']))
            else:
                columns['DNA Sequence count'].append(0)
            if 'RNA Sequence' in participants:
                columns['RNA Sequence count'].append(len(participants['RNA Sequence']))
            else:
                columns['RNA Sequence count'].append(0)
        else:
            # top level events don't have real event ids
            
            columns['Protein count'].append(0)
            columns['DNA Sequence count'].append(0)
            columns['RNA Sequence count'].append(0)
            
        columns['Pathway'].append(stId)
        trees.append(tree)
        names.append(event.name)
        columns['Reaction count'].append(event.reaction_count_deep())

    for i in range(0, len(trees)):
        label = trees[i] + names[i]
        columns['tree'].append(label + ' ' * (max_label - len(label)))        

    df = pd.DataFrame.from_dict(data=columns)

    if ctx.obj.output_directory:
        df.to_csv(
            path_or_buf=os.path.join(
                ctx.obj.output_directory, ctx.obj.file_prefix + 'pathway_nodes_with_reactions.csv'),
            mode='w'
        )

    if ctx.obj.show:
        print(df)

@cli.command()
@click.pass_context
def db_statistics(ctx):
    events = eventtree_(ctx)
    if events is None:
        print(f'No events at {ctx.obj.top_level}.')
        return
    orthologs = allOrthologs_(ctx)
    genes_by_reaction = {}
    for parent_name, parent_id, reaction_id, uniprot_id, rap_id, species_genes in orthologs:
        k = str(parent_id) + '-' + str(reaction_id)
        if k not in genes_by_reaction:
            genes_by_reaction[k] = set()
        genes_by_reaction[k] |= {uniprot_id}
    for name, tree in events.children.items():
        genes = set()
        for event, reaction in tree.all_reactions():
            k = str(event.stId) + '-' + str(reaction.stId)
            if k in genes_by_reaction:
                genes |= genes_by_reaction[k]
        print(name, tree.statistics(), len(genes))

if __name__ == '__main__':
    cli()
