# -*- coding: utf-8 -*-
"""

"""

import os
import sys
import zipfile

def check_file(zip_temp, tsa_str, nofile_list, exists_list):

    zip_path = zip_temp % tsa_str
    
    # Check if the zip file exsist
    if not os.path.exists(zip_path):
        #print 'No zip file for found for TSA %s with path %s' % (tsa_str, zip_path)
        # If not, append the TSA string to the list of invalid paths
        nofile_list.append(tsa_str)
        return None
    
    # Check it there is already a directory of the same name
    this_dir, this_file = os.path.split(zip_path)
    dirs = next(os.walk(this_dir))[1] # Returns just the dirs in this_dir
    bn = this_file.replace('.zip','')
    if bn in dirs:
        exists_list.append(tsa_str) #Append to list of existing dirs
        return None
    
    # If it hasn't already been zipped, return the path
    else:
        return zip_path
   

def unzip_file(path, out_dir):
    
    try:
        print 'Extracting %s to output directory %s' % (path, out_dir)
        z = zipfile.ZipFile(path)
        z.extractall(out_dir)
        z.close()
        
    except:
        print 'Problem extracting for ' + path
    
    

def main(tsa_txt, zip_template):
    
    df = pd.read_csv(tsa_txt, sep='\t', usecols=['tsa_str'], dtype=object)
    #zip_paths = [zip_template % tsa for tsa in df_tsa['tsa_str']]

    nofile_list = []
    exists_list = []
    z_paths = [check_file(zip_template, tsa, nofile_list, exists_list) for tsa in df['tsa_str']]
    z_paths = [z for z in z_paths if z]
    
    return z_paths

tsa_txt = '/vol/v2/stem/scripts/tsa_orwaca.txt'
zip_template = '/vol/v1/scenes/%s/outputs/nbr/history_m1o0.zip'

zips = main(tsa_txt, zip_template)
    
    
    