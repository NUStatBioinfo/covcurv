import argparse
import os
import pytest

def parse_args():
    """
    Obtain degnorm CLI parameters.

    :return: argparse.ArgumentParser object with runtime parameters required to run DegNorm pipeline.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--keep-output'
                        , action='store_true'
                        , required=False
                        , help='If specified keep covcurv_test output directory')

    args = parser.parse_args()
    return args


def main():

    args = parse_args()

    # determine whether or not to remove degnorm_test output directory after testing pipeline.
    os.environ['COVCURV_TEST_CLEANUP'] = 'True'
    if args.keep_output:
        os.environ['COVCURV_TEST_CLEANUP'] = 'False'

    # run existing unit tests.
    print('RUNNING covcurv TESTS...')
    tests_dir = os.path.dirname(os.path.abspath(__file__))
    pytest.main(['-x', tests_dir])


if __name__ == '__main__':
    main()