import multiprocessing as mp
import logging
import warnings
import sys
import re
import numpy as np
import subprocess
import os
from datetime import datetime
import time
import argparse
import pkg_resources
import gc


def configure_logger(output_dir, mpi=False):
    """
    Configure logger. Save log to file in output directory and route to stdout.

    :param output_dir: str path to output dir where covcurv.log file to be written.
    :param mpi: Bool is user running covcurv_mpi?
    """
    handlers = [logging.StreamHandler()]
    if output_dir:
        handlers += [logging.FileHandler(os.path.join(output_dir, 'covcurv.log'))]

    fmt = 'COVCURV (%(asctime)s) ---- %(message)s'
    if mpi:
        fmt = 'covcurv MPI (%(asctime)s) ---- %(message)s'

    logging.basicConfig(level=logging.DEBUG
                        , format=fmt
                        , handlers=handlers
                        , datefmt='%m/%d/%Y %I:%M:%S')


def welcome():
    """
    Welcome our user with DegNorm ascii art.
    """
    resources_dir = pkg_resources.resource_filename(__name__, 'resources')
    with open(os.path.join(resources_dir, 'welcome.txt'), 'r') as f:
        welcome = f.readlines()
        welcome += '\n' + 'version {0}'.format(pkg_resources.get_distribution('covcurv').version)

    logging.info('\n' + ''.join(welcome) + '\n'*4)


def create_output_dir(user_input=None):
    """
    Handle supplied output directory name for covcurv pipeline run.
    If no directory is supplied, use directory ./covcurv_<mmddYY>_<HHMMSS>
    If directory supplied is already a directory, use directory <supplied directory>/covcurv_<mmddYY>_<HHMMSS>
    If directory supplied is the name of a desired directory that does not exist already, use that.

    :param user_input: str user-supplied output directory
    :return: str decided upon and existent output directory
    """
    # if output directory is not specified, use default covcurv output directory naming convention.
    if not user_input:
        output_dir = os.path.join(os.getcwd(), 'covcurv_' + datetime.now().strftime('%m%d%Y_%H%M%S'))

    # if output directory is name of an existent directory, use default naming convention
    # and place output dir in the specified location.
    elif os.path.exists(user_input):
        output_dir = os.path.join(user_input, 'covcurv_' + datetime.now().strftime('%m%d%Y_%H%M%S'))

    # case when user has supplied a desired name for an output directory, needs to be created.
    else:
        output_dir = user_input

        # if just name of an output directory is given, place output dir in cwd.
        if os.path.dirname(user_input) == '':
            output_dir = os.path.join(os.getcwd(), user_input)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    return output_dir


def subset_to_chrom(df, chrom, reindex=False):
    """
    Subset a pandas.DataFrame with a 'chr' (chromosome) column to a particular
    chromosome. Reset the index if desired.

    :param df: pandas.DataFrame with a 'chr' column
    :param chrom: str or list; chromosome(s) with which to subset df
    :param reindex: bool indicator reset index? Default False: do not reset index.
    :return: pandas.DataFrame
    """
    if not isinstance(chrom, list):
        chrom = [chrom]

    if not reindex:
        sub_df = df[df['chr'].isin(chrom)]
    else:
        sub_df = df[df['chr'].isin(chrom)].reset_index(drop=True)

    if sub_df.empty:
        raise ValueError('Chromosome subsetting resulted in an empty DataFrame!')

    return sub_df


def max_cpu():
    """
    :return: int number of CPUs available on a node minus 1
    """
    return mp.cpu_count() - 1


def flatten_2d(lst2d, arr=True):
    """
    Flatten a 2-dimensional list of lists or list of numpy arrays into a single list or numpy array.

    :param lst2d: 2-dimensional list of lists or list of 1-d numpy arrays
    :param arr: Bool return numpy array or list?
    :return: 1-dimensional list or numpy array
    """
    lst1d = [elt for lst1d in lst2d for elt in lst1d]
    return np.array(lst1d) if arr else lst1d


def find_software(software='samtools'):
    """
    Determine if a software is in a $PATH.

    :return: True if software is in path, otherwise False
    """
    out = subprocess.run(['which {0}'.format(software)]
                         , shell=True)
    if out.returncode != 0:
        return False

    return True


def bai_from_bam_file(bam_file):
    """
    Simple helper function to change the file extension of a .bam file to .bai.
    """
    if not bam_file.endswith('.bam'):
        raise ValueError('{0} must have a .bam extension.'.format(bam_file))

    return bam_file[:-3] + 'bai'


def create_index_file(bam_file):
    """
    Create a BAM index file with samtools.

    :param bam_file: str realpath to .bam file for which we desire a .bai index file
    :return: str realpath to the created .bai file
    """
    bai_file = bai_from_bam_file(bam_file)

    # check if samtools is available in $PATH.
    samtools_avail = find_software('samtools')

    # Note that samtools is only available for Linux and Mac OS:
    # https://github.com/samtools/samtools/blob/develop/INSTALL
    if not samtools_avail:
        raise EnvironmentError('samtools not in found in PATH. samtools is required to convert {0} -> {1}'
                               .format(bam_file, bai_file))

    # run samtools index
    cmd = 'samtools index {0} {1}'.format(bam_file, bai_file)
    out = subprocess.run([cmd]
                         , shell=True)

    if out.returncode != 0:
        raise ValueError('{0} was not successfully converted into a .bai file'.format(bam_file))


def split_into_chunks(x, n):
    """
    Split a list into a set of n approximately evenly-sized sublists.

    :param x: iterable - a list, numpy array, tuple, etc. Not a generator.
    :param n: int number of chunks. If n >= len(x), split x into sublists of size 1.
    :return: list of lists
    """
    csize = int(np.ceil(len(x) / n))
    out = list()
    
    i = 0
    while i*csize < len(x):
        out.append(x[(i * csize):(i * csize + csize)])
        i += 1

    return out


def argparser():
    """
    Obtain covcurv CLI parameters.

    :return: argparse.ArgumentParser object with runtime parameters required to run coverage curve pipeline.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--bam-files'
                        , nargs='+'
                        , default=None
                        , required=False
                        , help='Aligned read files (can be multiple; separate with space) in .bam format. '
                               '.bam files may be for paired or single-end read experiments.')
    parser.add_argument('--bai-files'
                        , nargs='+'
                        , default=None
                        , required=False
                        , help='Input .bam index files (Can be multiple; separate with space). '
                               'Only use with --bam-files.'
                               '.bai files for files passed to --bam-files. Assumed .bai file order corresponds to '
                               'supplied .bam files. If --bam-files is supplied and not --bai-files, it is assumed '
                               'that you have .bai files in the same location as --bam-files, and with the same '
                               'basename, e.g. (A01.bam, A01.bai), (sample_27.bam, sample_27.bai)')
    parser.add_argument('--bam-dir'
                        , default=None
                        , required=False
                        , help='Input .bam/.bai data directory. Use instead of, or in addition to, specifying individual '
                               '.bam files. All .bam files (and .bai files with the same basename) in this directory '
                               'will be considered covcurv input.')
    parser.add_argument('-g'
                        , '--genome-annotation'
                        , type=str
                        , default=None
                        , required=False
                        , help='Genome annotation file.'
                               'Must have extension .gtf or .gff.'
                               'All non-exon regions will be removed, along with exons that appear in '
                               'multiple chromosomes and exons that overlap with multiple genes.')
    parser.add_argument('-o'
                        , '--output-dir'
                        , type=str
                        , default=None
                        , required=False
                        , help='Name for a covcurv output directory.'
                               'Name an directory for storing coverage data output by covcurv. '
                               'If not specified, directory ./covcurv_[mmddyy]_[HHMMSS] will be created.')
    parser.add_argument('-u'
                        , '--unique-alignments'
                        , action='store_true'
                        , help='Only retain reads that were uniquely aligned. All reads with '
                               'the flag "NH:i:<x>" with x > 1 will be dropped.')
    parser.add_argument('-p'
                        , '--proc-per-node'
                        , type=int
                        , required=False
                        , default=max_cpu()
                        , help='Number of processes to spawn per node, for within-node parallelization.'
                               'Defaults to the number of available cores (on the worker node) - 1. '
                               'If too high for a given node, reduce to the node\'s default.')
    parser.add_argument('-v'
                        , '--version'
                        , action='version'
                        , version='covcurv version {0}'.format(pkg_resources.get_distribution('covcurv').version)
                        , help='Display covcurv package version and exit.')
    parser.add_argument('-h'
                        , '--help'
                        , action='help'
                        , default=argparse.SUPPRESS
                        , help='covcurv: Per-gene read coverage matrix pipeline package.')

    return parser


def parse_args():
    """
    Parse command line arguments.

    :return: parsed argparse.ArgumentParser
    """
    parser = argparser()
    args = parser.parse_args()

    # --------------------------------------------------------- #
    # Parse and check runtime parameters.
    # --------------------------------------------------------- #

    # check validity of cores selection.
    max_ppn = max_cpu() + 1
    if args.proc_per_node > max_ppn:
        warnings.warn('{0} is greater than the number of available cores ({1}). Reducing to {2}'
                      .format(args.proc_per_node, max_ppn, max_ppn - 1))
        args.proc_per_node = max_ppn - 1

    # ensure that user has supplied fresh .bam/.bai files or a directory containing alignment files.
    if not args.bam_files and not args.bam_dir:
        raise ValueError('Must specify either --bam-files or --bam-dir.')

    if args.bam_files and args.bam_dir:
        raise ValueError('Both --bam-files and --bam-dir were provided! Not sure which data set to use.')

    # --------------------------------------------------------- #
    # Gather input RNA-Seq + genome annotation files.
    # --------------------------------------------------------- #

    # check validity of gene annotation file selection.
    if not args.genome_annotation:
        raise ValueError('If warm-start directory not specified, gene annotation file must be specified!')

    else:
        if not os.path.isfile(args.genome_annotation):
            raise FileNotFoundError('Gene annotation file {0} not found.'.format(args.genome_annotation))

    # check validity of file i/o selection.
    bam_files = list()
    bai_files = list()
    create_bai_files = list()

    # INPUT OPTION 1: a --bam-dir was specified.
    if args.bam_dir:

        # if user used both --bam-dir and --bam-files and/or --bai-files, yell at them. (only use one method).
        if args.bam_files is not None or args.bai_files is not None:
            raise ValueError('Do not specify both a --bam-dir and either --bam-files and/or --bai-files.'
                             'Use one input selection method or the other.')

        # check that the dir actually exists.
        if not os.path.isdir(args.bam_dir):
            raise NotADirectoryError('Cannot find --bam-dir {0}'.format(args.bam_dir))

        # scan directory for .bam files.
        for f in os.listdir(args.bam_dir):
            if f.endswith('.bam'):
                bam_files.append(os.path.join(args.bam_dir, f))

        # search for .bai files in the --bam-dir. If they don't exist, try to make them.
        for bam_file in bam_files:
            bai_file = re.sub('.bam$', '.bai', bam_file)

            # if .bai file under same basename as .bam file doesn't exist,
            # add it to list of .bai files that need to be created.
            if not os.path.isfile(bai_file):
                bai_files.append(bai_from_bam_file(bam_file))
                create_bai_files.append(bam_file)
            else:
                bai_files.append(bai_file)

    # INPUT OPTION 2: --bam-files and possibly --bai-files were specified.
    else:
        # ensure .bam files are actually .bam files.
        for bam_file in args.bam_files:
            if not bam_file.endswith('.bam'):
                raise ValueError('{0} is not a .bam file.'.format(bam_file))
            elif not os.path.isfile(bam_file):
                raise FileNotFoundError('Count not find .bam file {0}'.format(bam_file))
            else:
                bam_files.append(bam_file)

        # case where user has specified .bai files to accompany .bam files.
        if args.bai_files is not None:
            # if user has supplied an incorrect number of bai files, fail out.
            if len(args.bai_files) != len(bam_files):
                raise ValueError('Number of supplied .bai files does not match number of supplied .bam files.')

            # ensure .bai files are actually .bai files.
            for bai_file in args.bai_files:
                if not bai_file.endswith('.bai'):
                    raise ValueError('{0} is not a .bai file.'.format(bai_file))
                elif not os.path.isfile(bai_file):
                    raise FileNotFoundError('Count not find .bai file {0}'.format(bai_file))
                else:
                    bai_files.append(bai_file)

        # if user has not supplied any bai files: look for them under the same name
        # as each of the .bam files, or create new .bai files with samtools (if possible).
        else:
            for bam_file in bam_files:
                bai_file = re.sub('.bam$', '.bai', bam_file)

                # if .bai file under same name as .bam file doesn't exist,
                # add it to list of .bam files for which we need to create a .bai file.
                if not os.path.isfile(bai_file):
                    bai_files.append(bai_from_bam_file(bam_file))
                    create_bai_files.append(bam_file)
                else:
                    bai_files.append(bai_file)

    # ensure that input files are uniquely named.
    if len(bam_files) != len(set(bam_files)):
        raise ValueError('Supplied .bam files are not uniquely named!')

    # create parser attributes for bam/index files.
    args.bam_files = bam_files
    args.bai_files = bai_files
    args.create_bai_files = create_bai_files

    return args