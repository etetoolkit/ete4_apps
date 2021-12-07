#! /usr/bin/env python
"""
"""
import sys
import os
import datetime

from time           import time
from random         import random, choices
from hashlib        import md5

import h5py
import numpy as np

from ete4           import Tree


def printime(msg):
    print(msg +
          (' ' * (79 - len(msg.replace('\n', '')))) +
          '[' +
          str(datetime.datetime.fromtimestamp(time()).strftime('%Y-%m-%d %H:%M:%S')) +
          ']')



def seq2num_aa(seq):
    """
    Converts amino-acid symbols to 8-bits intergers.
    Consensus through bitwise AND will give rise to some groups based
    on BLOSUM62 matrix.
    """
    conv = {
        int('00110011', 2): 'C',  # \
        int('01100011', 2): 'S',  #  | - small and polar
        int('11000011', 2): 'R',  # /
        int('00110110', 2): 'P',  # \
        int('01100110', 2): 'A',  #  | - small and non-polar
        int('11000110', 2): 'G',  # /
        int('00111100', 2): 'N',  # \
        int('01101100', 2): 'D',  #  | _ polar or acidic
        int('11001100', 2): 'E',  #  |
        int('10011100', 2): 'Q',  # /
        int('00110101', 2): 'H',  # \
        int('01100101', 2): 'R',  #  | - basic
        int('11000101', 2): 'K',  # /
        int('00111001', 2): 'M',  # \
        int('01101001', 2): 'L',  #  |_ large and hydrophobic
        int('11001001', 2): 'I',  #  |
        int('10011001', 2): 'V',  # /
        int('00111010', 2): 'F',  # \
        int('01101010', 2): 'Y',  #  | - aromatic
        int('11001010', 2): 'W',  # /
        int('00000000', 2): 'X',
        int('11111111', 2): '-',
    }
    return np.array([conv[n] for n in seq], dtype='int8')


def _print_sequence_aa(sequence):
    rev_conv = {
        int('00110011', 2): 'C',  # \
        int('01100011', 2): 'S',  #  | - small and polar
        int('11000011', 2): 'R',  # /

        int('00110110', 2): 'P',  # \
        int('01100110', 2): 'A',  #  | - small and non-polar
        int('11000110', 2): 'G',  # /

        int('00111100', 2): 'N',  # \
        int('01101100', 2): 'D',  #  | _ polar or acidic
        int('11001100', 2): 'E',  #  |
        int('10011100', 2): 'Q',  # /

        int('00110101', 2): 'H',  # \
        int('01100101', 2): 'R',  #  | - basic
        int('11000101', 2): 'K',  # /

        int('00111001', 2): 'M',  # \
        int('01101001', 2): 'L',  #  |_ large and hydrophobic
        int('11001001', 2): 'I',  #  |
        int('10011001', 2): 'V',  # /

        int('00111010', 2): 'F',  # \
        int('01101010', 2): 'Y',  #  | - aromatic
        int('11001010', 2): 'W',  # /

        int('00000000', 2): 'X',
        int('11111111', 2): '-',

        # official conventions:
        int('01001100', 2): '-',  # negatively charged
        int('00101100', 2): 'B',  # N or D
        int('10001100', 2): 'Z',  # E or Q
        int('01001001', 2): 'J',  # I or L

        int('00100011', 2): 'p',  # \
        int('01000011', 2): 'p',  #  | - small and polar
        int('00000011', 2): 'p',  # /

        int('00100110', 2): 's',  # \
        int('01000110', 2): 's',  #  | - small and non-polar
        int('00000110', 2): 's',  # /

        int('00011100', 2): 'a',  # \
        # int('00101100', 2): 'a',  #  | 
        # int('01001100', 2): 'a',  #  | - polar or acidic
        # int('10001100', 2): 'a',  #  |
        int('00001100', 2): 'a',  # /

        int('00100101', 2): '+',  # \
        int('01000101', 2): '+',  #  | - basic
        int('00000101', 2): '+',  # /

        int('00011001', 2): 'l',  # \
        int('00101001', 2): 'l',  #  |
        # int('01001001', 2): 'l',  #  | - large and hydrophobic
        int('10001001', 2): 'l',  #  |
        int('00001001', 2): 'l',  # /

        int('00100101', 2): 'o',  # \
        int('01000101', 2): 'o',  #  | - aromatic
        int('00000101', 2): 'o',  # /

        # everything else is also X...
    }
    return ''.join(rev_conv.get(s, 'X') for s in sequence)


def seq2num_nt(seq):
    """
    Converts nucleotide symbols to 8-bits intergers
    """
    conv = {
        'A': int('10010001', 2),
        'C': int('10100010', 2),
        'G': int('01000111', 2),
        'T': int('01111000', 2),
        'W': int('00010000', 2),  # A or T
        'S': int('00000010', 2),  # C or G
        'R': int('00000000', 2),  # A or G
        'Y': int('00100000', 2),  # C or T
        'K': int('01000000', 2),  # G or T
        'M': int('10000000', 2),  # A or C
        'N': int('00000000', 2),
        '-': int('11111111', 2),
    }
    return np.array([conv[n] for n in seq], dtype='int8')


def consensus(sequences):
    """
    returns consensus sequence of an alignment as an array of int-8

    :param sequences: numpy 2D array of dtype int-8 (rows are individual 
       sequences)
    """
    return np.bitwise_and.reduce(sequences)


def _print_sequence_nt(sequence):
    rev_conv = {
        int('10010001', 2): 'A',
        int('10100010', 2): 'C',
        int('01000111', 2): 'G',
        int('01111000', 2): 'T',
        int('00010000', 2): 'W',  # A or T
        int('00000010', 2): 'S',  # C or G
        int('00000000', 2): 'R',  # A or G
        int('00100000', 2): 'Y',  # C or T
        int('01000000', 2): 'K',  # G or T
        int('10000000', 2): 'M',  # A or C
        int('00000000', 2): 'N',
        int('11111111', 2): '-',
    }
    return ''.join(rev_conv[s] for s in sequence)



def generate_random_alignment(tree, seq_len, fname, overwrite=False):
    """
    seq_len = 4_000_000_000 # genome size
    seq_len = 400_000_000   # chromosome size
    seq_len = 4_000_000     # big chunk
    seq_len = 40_000        # protein
    """

    if os.path.exists(fname) and not overwrite:
        return
    out = open(fname, 'w')
    for leaf_name in t.iter_leaf_names():
        seq = "".join("ATCG"[int(random() * 4)] for _ in range(seq_len))
        out.write(f'>{leaf_name}\n{seq}\n')
    out.close()


def fasta_reader(fasta, header_delimiter='\t'):
    with open(fasta) as fh:
        header = next(fh)[1:].strip().split(header_delimiter)
        seq = ''
        for l in fh:
            if l.startswith('>'):
                yield header, seq
                header = l[1:].strip().split(header_delimiter)
                seq = ''
                continue
            seq += l.strip()
        yield header, seq


def dump_h5tree(tree, h5out, fasta=None, overwrite=False, chunk_size=(1000, 100_000)):
    if os.path.exists(h5out) and not overwrite:
        return

    tree_len = len(tree)

    # with h5py.File('foo2.hdf5','w') as treef:
    treef = h5py.File(h5out, 'w', libver='latest') 
    #####################################
    # create GROUP for generic leaf data
    # leaf datasets can be represented globally, but each leaf represents 
    # a subset of it (i.e. an alignment)
    lg = treef.create_group("leaf_data")
    # - as all datasets here contain one entry per leaf, the attributes
    #   of this group will allow to to map entry numbers to the leaves
    leaf_number = 0
    for leaf_name in tree.iter_leaf_names():
        # check for leaf_name clash
        if leaf_name in lg.attrs:
            raise Exception(f'ERROR: Leaf name {leaf_name} found multiple times.')
        lg.attrs[leaf_name] = leaf_number
        leaf_number += 1
    # create a SEQUENCE dataset to store alignment
    if fasta:
        # get sequence length
        # initialize alignment dataset
        printime(' - Dumping FASTA to hdf5 alignment group')
        fr = fasta_reader(fasta)
        header, seq = next(fr)
        alignment = lg.create_dataset(
            'alignment', (tree_len, len(seq)), dtype='int8', 
            chunks=(min(tree_len, chunk_size[0]), min(len(seq), chunk_size[1]))) #, compression='gzip', compression_opts=7)
        alignment[lg.attrs[header[0]]] = seq2num_nt(seq)
        for header, seq in fr:
            alignment[lg.attrs[header[0]]] = seq2num_nt(seq)
    
    t1 = time()
    printime(' - Creating h5Tree groups')
    traverser = tree.traverse()
    h5root = next(traverser)
    root = treef.create_group("tree")
    h5root.add_prop('h5node', root)
    # store name, dist, support as h5py attributes
    for k, v in tree.props.items():
        root.attrs[k] = v
    # link to sequences
    #     root.attrs['sequence'] = alignment.ref  # store a reference to the alignment

    for node in traverser:
        if not node.is_leaf():  # we create a fancy hashed name
            hashed_name = md5(':'.join(node.get_leaf_names()).encode())
            g = node.up.props['h5node'].create_group(hashed_name.hexdigest())
        else:
            # group name is leaf name
            g = node.up.props['h5node'].create_group(node.name)
            # associate each leaf to each leaf-dataset (e.g. sequences)
            if treef['leaf_data'].keys():
                for ld in treef['leaf_data']:
                    g.attrs[ld] = lg.attrs[node.name]
        node.add_prop('h5node', g)
    return t1


def _load_tree_from_h5(tree, h5tree):
    """
    populate ete4 Tree using hdf5 group structure

    :param tree: ete4 Tree object (empty)
    :param h5tree: open instance of hdf5 file with tree
    """
    for h5n in h5tree:
        n = tree.add_child()
        n.name = h5n
        for k, v in h5tree[h5n].attrs.items():
            n.add_prop(k, v)
        n.add_prop('h5node', h5tree[h5n])
        _load_tree_from_h5(n, h5tree[h5n])


def load_h5tree(fname):
    """
    Load ete4 Tree from hdf5 attaching a h5node property 
    to each node.

    :returns: Tree object
    """
    h5node = h5py.File(fname)
    # create empty Tree to be populated from hdf5 file structure
    tree = Tree()
    tree.add_prop('h5node', h5node['tree'])
    _load_tree_from_h5(tree, h5node['tree'])
    return tree


times = {}

t = Tree()

tree_len = int(sys.argv[1])  # wanted random tree length
seq_len  = int(sys.argv[2])  # wanted random sequences length
chunk_size = int(sys.argv[3]), int(sys.argv[4])  # for storing alignment

t.populate(tree_len, names_library=map(str, range(tree_len + 1)))
base_name_fasta = f'{len(t)}leaves_{int(seq_len / 1_000)}Kbp'
base_name = f'{len(t)}leaves_{int(seq_len / 1_000)}Kbp_chunks{chunk_size[0]}-{chunk_size[1]}'
fasta  = f'alignment_{base_name_fasta}.fasta'

printime('Generate random alignment')
generate_random_alignment(t, seq_len, fasta)

printime('Create h5tree')
t0 = time()
t1 = dump_h5tree(t, f'h5tree_{base_name}.hdf5', fasta)
if t1 is None:
    t1 = time()
times['store hdf5 alignment'] = t1 - t0
times['create hdf5 tree'] = time() - t1
t1 = time()

del(t)

printime('Read h5tree')
h5tree = load_h5tree(f'h5tree_{base_name}.hdf5')
times['read hdf5 tree'] = time() - t1

printime('random access to h5tree sequences')
for _ in range(10):
    wanted_leaves = choices(h5tree.get_leaf_names(), k=1000)

    refs = sorted(set([(h5tree & l).props['h5node'].attrs['alignment']
                for l in wanted_leaves]))

    t1 = time()
    print(h5tree.props['h5node'].file['leaf_data']['alignment'][refs, 400:600])
    times.setdefault('random access alignment', []).append(time() - t1)

log = open(f'h5test_{base_name}.log', 'w')
for k in times:
    if k == 'random access alignment':
        log.write(f'{k}\t{",".join(str(t) for t in times[k])}\n')
    else:
        log.write(f'{k}\t{times[k]}\n')
log.close()

# clean 
# os.system(f'rm -f h5tree_{base_name}.hdf5')

printime('Done.')
