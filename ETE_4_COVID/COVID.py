

from ete4 import PhyloTree, SeqGroup

from ete4.smartview import TreeStyle, SeqMotifFace

from collections import OrderedDict

import argparse, sys, os

from importlib.machinery import SourceFileLoader


def check_arg(args = None):

    '''

    The function is used for parsing the input parameters form the command line

    using the standard python package argparse. The package itself is handling

    the validation and the return errors messages.

​

    '''

    parser = argparse.ArgumentParser(prog = 'multialignment.py',

                                    formatter_class = argparse.RawDescriptionHelpFormatter,

                                    description = 'Phylogenetic tree with Multialignment')

    parser.add_argument('-v','--version', action='version',version='%(prog)s 0.0.1')

    parser.add_argument('-t','--tree',required = True,

                                    help = 'Phylogenetic tree annotated')

    parser.add_argument('-m','--multialignment',required = True,

                                    help = 'Multialignment file with the reference genome')

    parser.add_argument('-c','--coordinates',required = True,

                                    help = 'Coordinates file')

    parser.add_argument('-l','--layouts',
                                    help= 'Layouts displayed in the graphical interface. The file MUST be named layouts_file and MUST be in the current directory')
    
    return parser.parse_args()
    


def check_if_file_exists(file_to_check):

    '''

    Description:

        The function will check if the file exists

    Input:

        file_to_check       # file name

    Return:

        True if exists , else False

    '''

    if os.path.isfile(file_to_check):

        return True

    else:

        return False
        
        

def read_coordinates_file (coord_file):

    dict_coord = OrderedDict ()

    with open (coord_file, 'r') as fh:

        for line in fh.readlines()[1:]:

            line=line.strip().split(',')

            if line == []:

                continue

            dict_coord[line[0]]=line[1:3]

    return dict_coord

            

def get_layout_alignment(reg, alg, col_name):

    def layout_fn(node):

        seq_format = 'seq'

        if node.is_leaf():

            seq = alg.get_seq(node.name)

            seqFace = SeqMotifFace(seq, seq_format=seq_format, padding_x=0, padding_y=0, width=len(seq)*15, height=20)

            node.add_face(seqFace, column = col_name, position = "aligned") 

        else:

            # I believe this gives you the first leaf

            first_leaf = next(node.iter_leaves())

            seq = alg.get_seq(first_leaf.name)

            seqFace = SeqMotifFace(seq, seq_format=seq_format, bgcolor="red", padding_x=0, padding_y=0, width=len(seq)*15, height=20)

            node.add_face(seqFace,  column = col_name, position = "aligned", collapsed_only=True) 

    layout_fn.__name__ = 'Alignment ' + reg

    layout_fn.contains_aligned_face = True

    return layout_fn
    

if __name__ == '__main__':    

    arguments = check_arg(sys.argv[1:])
    #print(arguments)

    if not check_if_file_exists(arguments.tree):

        print('Tree file does not exist\n')

        exit(2)

    if not check_if_file_exists(arguments.multialignment):

        print('Multialignment file does not exist\n')

        exit(2)

    if not check_if_file_exists(arguments.coordinates):

        print('Coordinates file does not exist\n')

        exit(2)


    from layouts_file import get_layout_clade, get_layout_lineage, get_layout_continent, get_layout_who_label,get_layout_nextstrain_clade

    #Load tree 

    t = PhyloTree(arguments.tree,format=1)

    #t.resolve_polytomy()

    ts = TreeStyle()

    

    #Load multialignment (Each item is a tuple with the sequence name, sequence, and sequence comments)

    alignments = SeqGroup(arguments.multialignment)



    #Get coordinates

    dict_coord = read_coordinates_file(arguments.coordinates)



    dict_region = OrderedDict ()

    for region in dict_coord.keys():

        dict_region[region] = SeqGroup()

        for alignment in alignments.iter_entries():

            # import pdb; pdb.set_trace()

            accession, seq = alignment[0], alignment[1][int(dict_coord[region][0]):int(dict_coord[region][1])]

            dict_region[region].set_seq(accession, seq)

    

    layouts = []

    col_name = 0

    for reg, alignment in dict_region.items():

        col_name+=1

        layouts.append(get_layout_alignment(reg, alignment,col_name))


    layouts=layouts+[get_layout_clade(), get_layout_lineage(), get_layout_continent(), get_layout_who_label(),get_layout_nextstrain_clade()]
    
    t.explore(tree_name="0", tree_style=ts, layouts=layouts)

