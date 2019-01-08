import pickle as pkl


def load_coverage_data(pkl_file):
    """
    Load a pickled up coverage data file (will be used by reticulate::source_python)

    :param pkl_file: str coverage data filename
    :return: loaded object stored in pkl_file
    """
    with open(pkl_file, 'rb') as f:
        cov_dat = pkl.load(f)

    return cov_dat