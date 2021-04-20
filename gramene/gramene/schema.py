import pandas as pd


def path_generator(path_str):
    length = len(path_str)
    a = 0
    # Skip any slashes at the start of the path.
    while a < length and path_str[a] == '/':
        a += 1
    b = a
    while b < length:
        if path_str[b] == '\\':
            # Skip escape sequences.
            b += 1
            if b < length and path_str[b] == '/':
                b += 1
        elif path_str[b] == '/':
            # Unescaped slashes.
            yield path_str[a:b]
            a = b + 1
            b = a
        else:
            b += 1
    if b > a + 1:
        yield path_str[a:b]


def normalize_path(path):
    if path is None:
        return ()
    if type(path) is tuple:
        return path
    if type(path) is list:
        return tuple(path)
    if type(path) is not str:
        raise ValueError(
            'Path must be a string, tuple, or list.  Received: ' + type(path).__name__)
    return tuple(path_generator(path))


class PathwayBase:
    def __init__(self, children):
        self.children = {}
        self.reactions = []
        self.black_box_events = []
        for item in children:
            if item['type'] in ('Pathway', 'TopLevelPathway'):
                pathway = Pathway(item)
                self.children[pathway.name] = pathway
            elif item['type'] == 'Reaction':
                self.reactions.append(Reaction(item))
            elif item['type'] == 'BlackBoxEvent':
                self.black_box_events.append(BlackBoxEvent(item))
            else:
                print(f'Unknown reaction type: {item["type"]}')

    def __str__(self):
        return self.name

    def __getitem__(self, path):
        path_tuple = normalize_path(path)
        if len(path_tuple) == 0:
            return self
        if path_tuple[0] in self.children:
            return self.children[path_tuple[0]][path_tuple[1:]]
        return None

    def keys(self):
        return self.children.keys()

    def walk(self, depth=0, prefix='', coming_siblings=0):
        ascii_tree = '+-- '
        children_prefix = prefix + '|   '
        if coming_siblings == 0:
            ascii_tree = '`-- '
            children_prefix = prefix + '    '
        yield self, depth, prefix + ascii_tree
        num_children = len(self.children)
        for item in self.children.values():
            num_children -= 1
            yield from item.walk(depth+1, children_prefix, num_children)

    def statistics(self):
        stats = {
            'reactions': len(self.reactions) + len(self.black_box_events),
            'pathways': 1,
        }
        reactions = 0
        pathways = 0
        for item in self.children.values():
            child_stats = item.statistics()
            stats['reactions'] += child_stats['reactions']
            stats['pathways'] += child_stats['pathways']
        return stats

    def to_data_frame(self):
        data = []
        for event, depth, tree, in self.walk():
            if event.name == 'TopLevel':
                continue
            data.append((event.stId, depth, tree, event.name,
                         event.species, event.diagram))
        return pd.DataFrame(
            data=data,
            columns=['stId', 'depth', 'tree', 'name', 'species', 'diagram']
        )

    def all_reactions(self):
        for item, depth, tree in self.walk():
            for event in item.reactions:
                yield item, event
            for event in item.black_box_events:
                yield item, event

    def all_reaction_ids(self):
        ids = set()
        for parent, reaction in self.all_reactions():
            ids.add(reaction.stId)
        return ids

    def reaction_count(self):
        return len(self.reactions) + len(self.black_box_events)

    def reaction_count_deep(self):
        count = 0
        for child in self.children.values():
            count += child.reaction_count_deep()
        return count + self.reaction_count()

    def pathway_nodes_with_reactions(self):
        for item, depth, tree in self.walk():
            if item.reaction_count() > 0:
                yield item


class EventsHierarchy(PathwayBase):
    def __init__(self, data):
        super().__init__(data)
        self.name = 'TopLevel'


class Pathway(PathwayBase):
    def __init__(self, data):
        super().__init__(data.get('children', []))
        self.stId = data['stId']
        self.name = data['name']
        self.species = data['species']
        self.diagram = data['diagram']


class Reaction:
    def __init__(self, data):
        self.stId = data['stId']
        self.name = data['name']
        if 'species' in data:
            self.species = data['species']
        else:
            self.species = ''

    def __str__(self):
        return self.name


class BlackBoxEvent(Reaction):
    def __init__(self, data):
        super().__init__(data)


class SimpleEntity:
    def __init__(self, id, data):
        self.id = id
