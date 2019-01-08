#!/usr/bin/env python3

# ---------------------------------------------------------------------------- #
# Launch covcurv web viz app
# ---------------------------------------------------------------------------- #

import os
import sys
import pkg_resources
import subprocess
import argparse


def check_covcurv_dir(data_dir):
    """
    Check that certain files exist within a DegNorm output directory.

    :param data_dir: str path to DegNorm pipeline run output directory
    :param file_names: str or list of str files output from DegNorm pipeline, contained in data_dir
    """
    if not os.path.isdir(data_dir):
        raise IOError('covcurv data directory {0} not found.'.format(data_dir))

    if not os.path.isfile(os.path.join(data_dir, 'gene_exon_metadata.csv')):
        raise IOError('Required covcurv output file "gene_exon_metadata.csv" is missing within {0}. '
                      'Ensure that {0} contains covcurv output data.'.format(data_dir))

def main():

    # ---------------------------------------------------------------------------- #
    # Parse (supplied) covcurv output directory
    # ---------------------------------------------------------------------------- #
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-d', '--data-dir'
                        , type=str
                        , default='.'
                        , required=False
                        , help='covcurv output directory file path')
    parser.add_argument('-p'
                        , '--python'
                        , type=str
                        , default=None
                        , help='Path to python3 binary. covcurv_app requires python3 in the background. '
                               'If not specified, it is assumed your default `python` command is a python3 interpreter.')
    parser.add_argument('-h'
                        , '--help'
                        , action='help'
                        , default=argparse.SUPPRESS
                        , help='covcurv_app: launch web application to visualize/export coverage matrices.')
    args = parser.parse_args()

    # run quality control checks.
    if not os.path.isdir(args.data_dir):
        raise IOError('covcurv data directory {0} not found.'.format(args.data_dir))

    req_files = [os.path.join(args.data_dir, x) for x in ['gene_exon_metadata.csv', 'read_counts.csv']]
    for req_file in req_files:
        if not os.path.isfile(req_file):
            raise IOError('Required covcurv output file {0} not found'.format(req_file))

    # ---------------------------------------------------------------------------- #
    # Welcome user.
    # ---------------------------------------------------------------------------- #
    # display covcurv logo
    sys.stdout.write('\nStarting covcurv app...\n')

    welcome_dir = pkg_resources.resource_filename(__name__, 'resources')
    with open(os.path.join(welcome_dir, 'welcome.txt'), 'r') as f:
        welcome = f.readlines()
        welcome += '\n' + 'version {0}'.format(pkg_resources.get_distribution('covcurv').version)

    sys.stdout.write('\n' + ''.join(welcome) + '\n'*2)

    sys.stdout.write('Data directory={0}...OK\n'.format(args.data_dir))
    sys.stdout.flush()

    # ---------------------------------------------------------------------------- #
    # Start process that launches shiny server (shiny::runApp),
    # route port info back to terminal
    # ---------------------------------------------------------------------------- #
    app_dir = pkg_resources.resource_filename(__name__, 'resources/shiny')
    cmd = 'Rscript -e "library(shiny); DATADIR <- \'{data_dir}\'; PYTHON <- \'{python}\'; runApp(\'{app_dir}\')"'\
        .format(data_dir=os.path.abspath(args.data_dir), python=args.python, app_dir=app_dir)

    subprocess.run([cmd]
                   , shell=True
                   , stdout=sys.stdout
                   , stderr=sys.stderr)


if __name__ == "__main__":
    main()