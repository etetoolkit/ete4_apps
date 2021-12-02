
from ete4 import Tree
import csv
import re
import time
import argparse, sys, os

start=time.time()

def check_arg(args = None):

    '''

    The function is used for parsing the input parameters form the command line

    using the standard python package argparse. The package itself is handling

    the validation and the return errors messages.

​

    '''

    parser = argparse.ArgumentParser(prog = 'props_tree.py',

                                    formatter_class = argparse.RawDescriptionHelpFormatter,

                                    description = 'Prepare metadata file needed for covid phylogeny ete4 view')

    parser.add_argument('-v','--version', action='version',version='%(prog)s 0.0.1')

    parser.add_argument('-i','--input',required = True,

                                    help = 'Tree input in newick format')

    parser.add_argument('-o','--output',required = True,

                                    help = 'Output tree with new properties added')

    parser.add_argument('-m','--metadata',required = True,

                                    help = 'Metadata use to annotate the tree (tsv file). It is a pre-configured metadata file: it is recommended to use prepare_metadata.py before using props_tree.py')
                                    
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


if __name__ == '__main__':

    arguments = check_arg(sys.argv[1:])

    if not check_if_file_exists(arguments.input):

        print('Tree file does not exist\n')

        exit(2)

    if not check_if_file_exists(arguments.metadata):

        print('Coordinates file does not exist\n')

        exit(2)


    #### 1.Load the tree #####
    print('Starting...')
    print('Loading tree...')
    t = Tree(arguments.input)

    leaf2metadata = {}

    counter=0
    print('ADDING NEW PROPERTIES...\n')
    print('Creating dictionary from metadata...\n')
    for line in open(arguments.metadata):
        counter+=1

        if counter==1:
            headers=line.strip().split("\t")
        else:

            #fields_dict={}
            fields = line.strip().split("\t")
            fields_dict={headers[i]: fields[i] for i in range(len(fields))}      
            leaf2metadata[fields_dict['Accession ID']] = fields_dict 

    #print(leaf2metadata)
    #--Open a file called exceptions to store the sequences that cannot be annotated because they have been removed from gisaid (at least from metadata)
    exceptions=open("./exceptions.tsv", 'a')

    print('Iterating trough the leaves...\n')    
    number_leaves=0
    for leaf in t.iter_leaves():
        number_leaves+=1
        #print(number_leaves)
        leaf_name=leaf.name #The tree must have the accession id as the leaf name
        #leaf_name=leaf.name.split('|')[1] #This is for extracting the accession id if the leaf.name is of the form "XXXXX|EPI_ISL_NNNNNN|XXXXXX..."
        
        if re.match('^EPI_ISL*',leaf_name):
            accession_id=leaf_name #leaf_name.split("|")[1]
            #print('Accession ID',accession_id)
        else:
            print('Accession ID format is not appropiated, it must be: EPI_ISL_NNNNNN')
            print('Aborting')
            
        
        if accession_id in leaf2metadata.keys():
            data=leaf2metadata[accession_id]
            #print(data)
            data = {key : value.replace(':','') for key,value in data.items()} #directly format input data just in case it contains any : or = symbols
            data = {key : value.replace('=','') for key,value in data.items()}
            
            for key in data.keys(): #Be carefull, because the name of the geolocation property MUST be latlong.
                    leaf.add_prop(key.replace(' ','_').lower(),data[key])
                    
            leaf.add_prop('country',leaf.props['location'].split('/')[1].strip())
        else:
            #print('No such accession id in metadata')
            exceptions.write(f"\n{accession_id}\n")
        #{node X: {'G':10, 'GK':40, 'L':80...}} 

    print('Iterating trough nodes adding properties...\n')
    number_nodes=0
    for node in t.traverse("postorder"):
        number_nodes+=1
        #print(number_nodes)
        node_dict={}
        regions={} 
        next_dict={}
        who_dict={}
        if node.is_leaf()==False:
            for l in node.iter_leaves():
                #print(l.props)
                if 'clade' in l.props:
                    if l.props['clade'] in node_dict.keys():
                        node_dict[l.props['clade']]+=1
                    else:
                        node_dict[l.props['clade']]=1
                if 'location' in l.props:
                    if l.props['location'].split('/')[0].strip() in regions.keys():
                        regions[l.props['location'].split('/')[0].strip()]+=1
                    else:
                        regions[l.props['location'].split('/')[0].strip()]=1
                if 'nextstrain_clade' in l.props:
                    if l.props['nextstrain_clade'] in next_dict.keys():
                        next_dict[l.props['nextstrain_clade']]+=1
                    else:
                        next_dict[l.props['nextstrain_clade']]=1
                if 'who_label' in l.props:
                    if l.props['who_label'] in who_dict.keys():
                        who_dict[l.props['who_label']]+=1
                    else:
                        who_dict[l.props['who_label']]=1



            node.add_prop('clade_frequency',[f'{key}_{value}' for key,value in node_dict.items()])
            node.add_prop('region_frequency',[f'{key}_{value}' for key,value in regions.items()])
            node.add_prop('nextstrain_frequency',[f'{key}_{value}' for key,value in next_dict.items()])
            node.add_prop('who_label_frequency',[f'{key}_{value}' for key,value in who_dict.items()])



    #--Write new tree with the properties added to a file:
    print('Writing new tree to file...')
    t.write(outfile=arguments.output, properties=[])
    print('\nIt took', round(time.time()-start,2), 'seconds.')