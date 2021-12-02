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


def seq2num(seq):
    """
    Converts nucleotide symbols to 8-bits intergers

    bits correspond to:
    A
    C
    G
    T
    N
    """
    conv = {
        'A': int('10000', 2),
        'C': int('00010', 2),
        'G': int('00110', 2),
        'T': int('11000', 2),
        '-': int('00000', 2),
    }
    return np.array([conv[n] for n in seq], dtype='int8')


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


def create_h5tree(tree, h5out, fasta=None, overwrite=False):
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
        fr = fasta_reader(fasta)
        header, seq = next(fr)
        alignment = lg.create_dataset('alignment', (tree_len, len(seq)), dtype='int8', 
                                      chunks=(tree_len, min(len(seq), 10_000))) #, compression='gzip', compression_opts=7)
        alignment[lg.attrs[header[0]]] = seq2num(seq)
        for header, seq in fr:
            alignment[lg.attrs[header[0]]] = seq2num(seq)

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


t = Tree()

tree_len = int(sys.argv[1])  # wanted random tree length
seq_len  = int(sys.argv[2])  # wanted random sequences length

t.populate(tree_len, names_library=map(str, range(tree_len + 1)))

fasta  = f'alignment_{len(t)}leaves_{int(seq_len / 1_000_000)}Mbp.fasta'

printime('Generate random alignment')
generate_random_alignment(t, seq_len, fasta)

printime('Create h5tree')
create_h5tree(t, f'h5tree_{len(t)}leaves_{int(seq_len / 1_000_000)}Mbp.hdf5', fasta)

printime('Read h5tree')
h5tree = load_h5tree(f'h5tree_{len(t)}leaves_{int(seq_len / 1_000_000)}Mbp.hdf5')

printime('random access to h5tree sequences')
wanted_leaves = choices(h5tree.get_leaf_names(), k=100)

refs = sorted(set([(h5tree & l).props['h5node'].attrs['alignment']
               for l in wanted_leaves]))

print(h5tree.props['h5node'].file['leaf_data']['alignment'][refs, 400:600])

printime('Done.')
