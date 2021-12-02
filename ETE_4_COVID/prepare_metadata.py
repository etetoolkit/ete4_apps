
from countryinfo import CountryInfo
from geopy.geocoders import Nominatim
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

    parser = argparse.ArgumentParser(prog = 'prepare_metadata.py',

                                    formatter_class = argparse.RawDescriptionHelpFormatter,

                                    description = 'Prepare metadata file needed for covid phylogeny ete4 view')

    parser.add_argument('-v','--version', action='version',version='%(prog)s 0.0.1')

    parser.add_argument('-i','--input',required = True,

                                    help = 'Metadata input file. Mandatory fields: Accession ID, Pango lineage, Location, Clade [from gisaid]')

    parser.add_argument('-o','--output',required = True,

                                    help = 'Metadata output file.')

    parser.add_argument('-t','--tsv',required = True,

                                    help = 'tsv correlation file (pango lineage, nextstrain clades, who labels), provided: corr_nomenclatures.py')
                                    
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

    
    geolocator = Nominatim(user_agent="carlaalvi7@gmail.com", timeout=10) #geolocator now with my email, think about changing it...
    arguments = check_arg(sys.argv[1:])

    if not check_if_file_exists(arguments.input):

        print('Tree file does not exist\n')

        exit(2)

    if not check_if_file_exists(arguments.tsv):

        print('Coordinates file does not exist\n')

        exit(2)

##--Reading corr_nomenclatures.tsv file as dictionary:
    print('Reading corr_nomenclatures.tsv file as dictionary...')
    nom_dict={}
    with open(arguments.tsv,'r') as f:
        for line in f:
            line=line.strip()
            (pango,nextstrain,who)=line.split('\t')
            nom_dict[pango]=[nextstrain,who]

    ##--Read input file and open output file using csv:
    print('Read input file and open output file using csv')
    with open(arguments.input,'r') as csvinput:
        with open(arguments.output, 'w') as csvoutput:
            writer = csv.writer(csvoutput, lineterminator='\n', delimiter='\t')
            reader = csv.reader(csvinput, delimiter='\t')
            counter=0
            all = []
            already_loc={} #e.g.: {'Spain': 'lat,long', 'China': 'lat,long' ...} lat,long de su CAPITAL
            
            ##--Iterate throuth input metadata rows to append the new information as new columns with headers: Latlong, Nextstrain clade and WHO label.
            print('Runing program...')
            for row in reader: 
                counter+=1
                #print(counter)

                ##--Write the headers
                if counter==1:
                    row.append('Latlong')
                    row.append('Nextstrain clade')
                    row.append('WHO label')
                    all.append(row)
                    location_pos=[idx for idx,r in enumerate(row) if r.casefold()=='Location'.casefold()][0]
                    pango_id_pos=[idx for idx,r in enumerate(row) if r.casefold()=='Pango lineage'.casefold()][0]

                ##--Write the new information for each sequence
                else:

                    ##--Search countries and its capitals to get the latitude and longitude of the lasts: 

                    country=row[location_pos].split("/")[1].strip() 
                    
                    if country not in already_loc.keys():
                        try:
                            capital=CountryInfo(country).capital()
                            loc = geolocator.geocode(capital+','+ country)                            

                        except Exception: #Making this exception because gambia was the unique country written in an "incorrect" way,
                                        #--so one has to be sure that the name of the countries is "standard english"
                            if country=='Gambia':
                                capital=CountryInfo('The Gambia').capital()
                                loc = geolocator.geocode(capital+','+ 'The Gambia')

                        #--Capital geolocation MUST be in the format lat,long
                        cap_geoloc=f"{loc.latitude},{loc.longitude}"
                        already_loc[country]=cap_geoloc

                    else:
                        cap_geoloc=already_loc[country]

                    row.append(cap_geoloc)
                    
                    ##--Write nextstrain clades and who labels: 
                    
                    pango_id=row[pango_id_pos].strip() 
                    #print('PANGO ID: ', pango_id)
                    #--If the pango_id that appears in the metadata is a subclade of the pango lineage in the tsv file,
                    #--it should still be the same nextstrain clade, because they are wider than pango lineages,
                    #--and the same happens with WHO labels.
                    pango_id=[pango_key for pango_key in nom_dict.keys() if re.match(f'^{re.escape(pango_key)}', pango_id)]
                        #--re.escape transforms the special regex characters like '.' to match them as string, and not as special,
                        #--for example, it transforms the dot '.' into '\.', that is used to specifically find dots, not any character.
                    
                    if pango_id!=[]:
                        try:
                            pango_id=pango_id[0]
                            nexstrain_clade=nom_dict[pango_id][0]
                            WHO_label=nom_dict[pango_id][1]

                            row.append(nexstrain_clade)
                            row.append(WHO_label)  

                        except Exception:

                            ##--All AY... pango ids corresponde to Delta variant, so make sure they are correctly annotated
                            if re.match('AY', pango_id):
                                nexstrain_clade=nom_dict['AY'][0]
                                WHO_label=nom_dict['AY'][1]
                                row.append(nexstrain_clade)
                                row.append(WHO_label)

                            ##--If there is not an actual correlation between pango and nextstrain clade or who label, annotate as uncorrelated
                            else:
                                nexstrain_clade='uncorrelated'
                                WHO_label='uncorrelated'
                                row.append(nexstrain_clade)
                                row.append(WHO_label)
                        
                        all.append(row)
                    else:
                        try:
                            pango_id=row[pango_id_pos].strip() 
                            nexstrain_clade=nom_dict[pango_id][0]
                            WHO_label=nom_dict[pango_id][1]

                            row.append(nexstrain_clade)
                            row.append(WHO_label)  

                        except Exception:

                            ##--All AY... pango ids corresponde to Delta variant, so make sure they are correctly annotated
                            if re.match('AY', pango_id):
                                nexstrain_clade=nom_dict['AY'][0]
                                WHO_label=nom_dict['AY'][1]
                                row.append(nexstrain_clade)
                                row.append(WHO_label)

                            ##--If there is not an actual correlation between pango and nextstrain clade or who label, annotate as uncorrelated
                            else:
                                nexstrain_clade='uncorrelated'
                                WHO_label='uncorrelated'
                                row.append(nexstrain_clade)
                                row.append(WHO_label)
                        
                        all.append(row)


            print('Writing info to output file...')
            writer.writerows(all)

print(time.time()-start, ' seconds')