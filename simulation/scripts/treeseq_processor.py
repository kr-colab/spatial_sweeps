import os, time, glob, argparse
import pandas as pd
import gzip
from spatial_sweeps.simulation.scripts.frequency_area import process_tree, tree_frequency_area

parser=argparse.ArgumentParser()
parser.add_argument('--infile')
parser.add_argument('--NW')
parser.add_argument('--SC')
parser.add_argument('--ID')
args=parser.parse_args()

def process_treeseq(filename, simpath):
    """
    Processing details
    """
    tree = process_tree(filename, f'{simpath}/simulation.log')
    df = tree_frequency_area(tree)
    df.to_csv(filename.replace('.trees','_frequency_area.txt'), sep='\t')
    #with open(filename, 'rb') as f_in, gzip.open(filename+'.gzip', 'wb') as f_out:
    #    f_out.writelines(f_in)
    #os.remove(filename)
    # open and process fil
        # delete file
        # save cool stuff

NW=args.NW
SC=args.SC
ID=args.ID
infile=args.infile
simpath=f'out/sim_{ID}_s_{SC}_Nw_{NW}'

process_treeseq(infile, simpath)
