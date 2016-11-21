# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 16:29:56 2016

@author: shooper
"""
import os
import fnmatch
import glob
import subprocess
import sys
import pandas as pd

def read_params(params):
    '''
    Read parameter file and parse into dictionary
    '''
    if not os.path.exists(params):
        print 'Parameter file given does not exist:\n%s' % params
        return None
    
    d = {}
    try:
        with open(params) as f:
            input_vars = [line.split(";") for line in f]
    except: 
        print 'Problem reading parameter file:\n%s' % params
        return None
    
    for var in input_vars:
        if len(var) == 2:
            d[var[0].replace(" ", "")] =\
            '"%s"' % var[1].strip(" ").replace("\n", "")

    print 'Parameters read from:\n%s\n' % params
    
    return d
    
    
    
    
def main(params):

    inputs = read_params(params)

    for i in inputs:
        exec ("{0} = str({1})").format(i, inputs[i])
    n_jobs = int(n_jobs)
    
    # Get the full path of all files in the basepath directory tree that match
    #   the search_str 
    matched = []
    for root, dirs, files in os.walk(basepath, followlinks=True):
        these_paths = [os.path.join(root, f) for f in files]
        these_paths = fnmatch.filter(these_paths, search_str)
        matched.extend(these_paths)
    
    df = pd.read_csv(tsa_txt, sep='\t')
    tsa_ls = ['0' + str(tsa) for tsa in df.tsa_id]
    #unmatched_tsas = [tsa for tsa in tsa_ls for f in matched if tsa not in f]
    #unmatched_tsas = [tsa for tsa in tsa_ls if tsa not in matched]
    # Get all tsas that are not in one of the matched files
    #unmatched_tsas = tsa_ls
    for f in matched:
        [tsa_ls.remove(tsa) for tsa in tsa_ls if tsa in f]# and tsa not in unmatched_tsas]
        #unmatched_tsas.extend(ls)
        
    these_tsas = tsa_ls[:n_jobs]
    for tsa in these_tsas:
        
        try:
            this_bp = os.path.join(basepath, '%s/outputs/nbr/' % tsa)
            if not os.path.exists(this_bp):
                print 'No scenes output directory found for tsa %s \nContinuing...\n' % tsa
                continue
            vertyrs_path = glob.glob(os.path.join(this_bp, '*_vertyrs.bsq'))[0]
            segmse_path  = glob.glob(os.path.join(this_bp, '*_segmse.bsq'))[0]
    
            # Get basename from segmse so you know which file the output was derived from
            bn = os.path.basename(segmse_path).replace('_segmse.bsq', '')
            out_dir = os.path.join(this_bp, 'derived_outputs')
            if not os.path.exists(out_dir): os.mkdir(out_dir)
            out_path = os.path.join(out_dir, '%s_mse_split.bsq' % bn)
        
        except Exception as e:
            print e
            continue
        
        #mse_split.py {vertyrs_path} {mse_path} {output_path}
        py_path = '/vol/v1/general_files/script_library/magnitude_by_year/mse_split.py'
        cmd = 'python {0} {1} {2} {3}'.format(py_path, vertyrs_path, segmse_path, out_path)
        print 'Processing scene %s\n' % tsa
        subprocess.call(cmd, shell=True)# this doesn't submit multiple processes at once. Consider using multiprocessing.Process: https://docs.python.org/2/library/multiprocessing.html#multiprocessing.Process

if __name__ == '__main__':
    params = sys.argv[1]
    sys.exit(main(params))

'''tsa_txt = '/vol/v2/stem/scripts/tsa_orwaca.txt'
basepath = '/vol/v1/scenes'
search_str = '*_mse_split.bsq'
n_jobs = 1'''
#params = '/vol/v2/stem/scripts/batch_split_mse_params.txt'
#main(params)
    