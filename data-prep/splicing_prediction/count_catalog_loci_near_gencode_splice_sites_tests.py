import unittest

import numpy as np

from count_catalog_loci_near_gencode_splice_sites import (
    compute_distances_to_nearest_splice_site, compute_intron_boundaries)


class ComputeIntronBoundariesTests(unittest.TestCase):

    def test_unsorted_exons(self):
        self.assertEqual(compute_intron_boundaries([(300, 400), (100, 200)]), [(201, 299)])

    def test_single_exon_has_no_introns(self):
        self.assertEqual(compute_intron_boundaries([(100, 200)]), [])

    def test_adjacent_exons_have_no_intron(self):
        self.assertEqual(compute_intron_boundaries([(100, 200), (201, 300)]), [])


class ComputeDistancesToNearestSpliceSiteTests(unittest.TestCase):

    def test_distances(self):
        sites = np.array([100, 200])
        # Tracts (1-based bases): 91-95, 96-100 (contains 100), 101-110, 150-151, 211-220
        start_0based = np.array([90, 95, 100, 149, 210])
        end = np.array([95, 100, 110, 151, 220])
        np.testing.assert_array_equal(
            compute_distances_to_nearest_splice_site(start_0based, end, sites), [5, 0, 1, 49, 11])

    def test_no_sites(self):
        np.testing.assert_array_equal(
            compute_distances_to_nearest_splice_site(np.array([10]), np.array([20]), np.array([], dtype=int)), [-1])


if __name__ == "__main__":
    unittest.main()
