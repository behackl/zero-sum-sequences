from importlib.metadata import version

import zero_sum_sequences


def test_package_version_matches_distribution_metadata():
    assert zero_sum_sequences.__version__ == version("zero-sum-sequences")
