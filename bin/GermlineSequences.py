#!/usr/bin/python3
"""
Adds germline sequences to tsv database file
Parses the imgt database reference files and creates new column for each germline_(v|d|j)_call specified in the tsv
Does not check if the all the reference files genes are present in the the germline_(v|d|j)_call columns already present in the tsv, to allow flexibility in the input for references (there can be more than needed)
"""

print("Debut")

# Info
__author__ = 'Chloe Taurel, Gael Millot'

# Imports
import os
from argparse import ArgumentParser, RawTextHelpFormatter
from collections import OrderedDict
from itertools import groupby
from textwrap import dedent
from time import time

import re
import pandas as pd

from changeo.IO import readGermlines




def addGermlineSequences(tsv_file, references, gaps, nogaps):
    """
    Write germline sequences for each vdj gene to tab-delimited database file

    Arguments:
      tsv_file : input tab-delimited database file.
      references : folders and/or files containing imgt repertoire data in FASTA format.
      gaps : if true, the new germline_._seq columns will contain gaps present in IMGT database
      nogaps : if true, new columns will not contain the gaps present in IMGT database, and will be called germline_._seq_no_gaps.

    Returns:
      out_file: names of the tsv output file with new columns.
    """
    dict_ref = readGermlines(references)

    df = pd.read_csv(tsv_file, sep='\t')
    
    # Get the input fields
    matching_cols = [col for col in df.columns if re.fullmatch(r"germline_._call", col)]
    
    if (len(matching_cols) < 1):
        raise ValueError(f"\n\n================\n\nERROR IN ", {__name__}, "\nTsv input file : ", {tsv_file}, " \ndoes not contain any column name matching pattern \"germline_._call\"\n\n================\n\n")

    for col in matching_cols:
        x = col[len("germline") : -len("_call")]
        seq_col = f"clonal_germline{x}_seq_with_gaps"
        seq_col_nogap = f"clonal_germline{x}_seq_no_gaps"

        unique_vals = df[col].dropna().unique()
        missing = [val for val in unique_vals if val not in dict_ref]
        if missing:
            raise ValueError(
                f"\n\n================\n\nERROR IN {__name__}\nTsv input file: {tsv_file}\n"
                f"The following value(s) in column '{col}' do not match any key in dict_ref:\n"
                f"{missing}\n\n================\n\n"
            )

        mapped_with_gaps = df[col].map(dict_ref)
        if gaps:
            df[seq_col] = mapped_with_gaps

        if nogaps:
            # Remove all '.' even if none exist (harmless)
            df[seq_col_nogap] = mapped_with_gaps.apply(lambda s: s.replace('.', '') if isinstance(s, str) else s)



    dir_name = os.path.dirname(tsv_file)
    base_name = os.path.basename(tsv_file)
    root, ext = os.path.splitext(base_name)
    out_file = os.path.join(dir_name, root + "_germ-seq" + ext)

    # --- save with all logical values in uppercase ---
    df = df.replace({True: "TRUE", False: "FALSE"})

    df.to_csv(out_file, sep='\t', index=False)
    return out_file






def getArgParser():
    """
    Defines the ArgumentParser

    Arguments:
    None

    Returns:
    an ArgumentParser object
    """
    # Define input and output field help message
    fields = dedent(
             '''
             output files:
                 germ-pass
                    tsv file with assigned germline sequences for each germline_(v|d|j)_call column found in the input file.

             required fields:
                 germline_v_call, germline_j_call or germline_d_call (at least one)

             output fields:
                 germline_v_seq, germline_d_seq, germline_j_seq (at least one)
              ''')
    # Define argument parser
    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter, epilog=fields)

    # Germlines arguments
    group = parser.add_argument_group('germline sequence finder arguments')
    group.add_argument("-i", '--input', action='store', dest='tsv_input_file', required=True,
                        help='''Tsv file containing at least the following fields: germline_v_call, 
                        germline_j_call, and germline_d_call (only if a D reference file is specified)''')
    group.add_argument('-r', '--ref', nargs='+', action='store', dest='references', required=True,
                        help='''List of imgt folders and/or fasta files (with .fasta, .fna or .fa
                         extension) with germline sequences. Ex : imgt_mouse_IGHJ.fasta''')
    group.add_argument('-g', '--gaps', action='store_true', help='''Create new 'germline_._seq' columns containing the germline sequences with alignment gaps ('.') as present in the IMGT database.
                         If both --gaps and --nogaps are specified, both types of columns are added.
                         If neither --gaps nor --nogaps is specified, both types of columns are added by default.''')
    group.add_argument('-n', '--nogaps', action='store_true', help='''Create new 'germline_._seq_no_gaps' columns containing the germline sequences with all alignment gaps ('.') removed.
                         If both --gaps and --nogaps are specified, both types of columns are added.
                         If neither --gaps nor --nogaps is specified, both types of columns are added by default.''')
    return parser



if __name__ == '__main__':
    """
    Parses command line arguments and calls main
    """

    # example for testing :
    
    # tsv_file = '/home/ctaurel/python_scripts/1_productive_seq_clone-pass_germ-pass.tsv'
    # references = ['/mnt/c/Users/ctaurel/Documents/imgt_ref_sequences_database/imgt/mouse/vdj/imgt_mouse_IGHJ.fasta', '/mnt/c/Users/ctaurel/Documents/imgt_ref_sequences_database/imgt/mouse/vdj/imgt_mouse_IGHV.fasta', '/mnt/c/Users/ctaurel/Documents/imgt_ref_sequences_database/imgt/mouse/vdj/imgt_mouse_IGHD.fasta']


    # Parse command line arguments
    parser = getArgParser()
    args = parser.parse_args()

    gaps = args.gaps
    nogaps = args.nogaps

    # Handle the default (neither argument provided)
    if not gaps and not nogaps:
        gaps = True
        nogaps = True


    # Check that reference files exist
    for f in args.references:
        if not os.path.exists(f):
            parser.error(f'Germline reference file or folder {f} does not exist.')

    if not os.path.exists(args.tsv_input_file):
        parser.error('Tsv file %s does not exist.' % tsv_file)


    addGermlineSequences(args.tsv_input_file, args.references, gaps, nogaps)
