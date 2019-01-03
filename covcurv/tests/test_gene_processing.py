import pytest
import os
from pandas import DataFrame
from covcurv.gene_processing import GeneAnnotationProcessor

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------- #
# GeneAnnotationProcessor tests
# ----------------------------------------------------- #

@pytest.fixture
def gtf_processor_setup():
    gtf_file = os.path.join(THIS_DIR, 'data', 'human_chr4_7.gtf')
    gtf_processor = GeneAnnotationProcessor(gtf_file
                                            , verbose=False)
    return gtf_processor


def test_gtf_processor_load(gtf_processor_setup):
    exons_df = gtf_processor_setup.load()
    assert isinstance(exons_df, DataFrame)
    assert not exons_df.empty


def test_gtf_processor_run(gtf_processor_setup):
    reqd_cols = ['chr', 'gene', 'gene_start', 'gene_end', 'start', 'end']
    exons_df = gtf_processor_setup.run()
    assert isinstance(exons_df, DataFrame)
    assert not exons_df.empty
    assert all([col in exons_df.columns.tolist() for col in reqd_cols])