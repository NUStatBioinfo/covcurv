import pytest
import os
import subprocess
import shutil
from random import choice

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------- #
# Whole-pipeline test
# ----------------------------------------------------- #
@pytest.fixture
def run_setup(request):

    # create temporary directory for test.
    salt = ''.join([choice(['a', 'b', 'c', 'd', 'e', 'f', 'g', '1', '2', '3', '4', '5', '6', '7']) for _ in range(20)])
    outdir = os.path.join('covcurv_test_' + salt)
    os.makedirs(outdir)

    def teardown():

        # determine whether or not to remove pipeline test output.
        if os.environ['COVCURV_TEST_CLEANUP'] == 'True':
            shutil.rmtree(outdir)

    request.addfinalizer(teardown)

    return outdir


def test_pipeline(run_setup):

    # run covcurv command with test data, with downsampling.
    cmd = 'covcurv --bam-files {0} {1} --bai-files {2} {3} -g {4} -o {5}'\
        .format(os.path.join(THIS_DIR, 'data', 'T2_FF_small.bam')
                , os.path.join(THIS_DIR, 'data', 'T2_FFPE_small.bam')
                , os.path.join(THIS_DIR, 'data', 'T2_FF_small.bai')
                , os.path.join(THIS_DIR, 'data', 'T2_FFPE_small.bai')
                , os.path.join(THIS_DIR, 'data', 'human_chr4_7.gtf')
                , run_setup)

    out = subprocess.run([cmd]
                         , shell=True
                         , stderr=subprocess.PIPE)
    assert out.returncode == 0

    # # run covcurv_mpi command with test data, without downsampling, without specifying .bai files either.
    # cmd = 'mpiexec -n 2 covcurv_mpi --bam-files {0} {1} -g {2} -o {3} -p 2 --plot-genes {4} --nmf-iter 50'\
    #     .format(os.path.join(THIS_DIR, 'data', 'hg_small_1.bam')
    #             , os.path.join(THIS_DIR, 'data', 'hg_small_2.bam')
    #             , os.path.join(THIS_DIR, 'data', 'chr1_small.gtf')
    #             , run_setup)

    # # if MPI available and > 1 node available, run pipeline test for covcurv_mpi.
    # try:
    #     from mpi4py import MPI
    #
    #     # try running covcurv_mpi test if there are > 1 nodes available. O/w skip it.
    #     try:
    #         out = subprocess.run([cmd]
    #                              , shell=True
    #                              , stderr=subprocess.PIPE)
    #
    #         assert out.returncode == 0
    #
    #     except RuntimeError:
    #         pass
    #
    # except ImportError:
    #     pass