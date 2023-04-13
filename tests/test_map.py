import unittest

import numpy as np

from geosample import MapSamples


class TestMap(unittest.TestCase):
    def test_map13(self):
        error_matrix = np.array(
            [[97, 0, 3], [3, 279, 18], [2, 1, 97]], dtype='float64'
        )

        sa = MapSamples(
            class_area=[22353, 1122543, 610228],
            error_matrix=error_matrix,
        )

        # Proportional class area
        self.assertTrue(
            np.allclose(sa.w, np.array([0.01273585, 0.63958045, 0.3476837]))
        )
        # Stratified unbiased area estimate
        self.assertTrue(
            np.allclose(sa.p_j, np.array([0.02570326, 0.59828666, 0.37601009]))
        )
        # Standard error
        self.assertTrue(
            np.allclose(
                sa.s_pj, np.array([0.00612572, 0.01005743, 0.01061797])
            )
        )
        # Standard error of the error-adjusted area conf
        self.assertTrue(
            np.allclose(
                sa.s_aj_conf,
                np.array([21072.36561, 34597.37001191, 36525.60632938]),
            )
        )
        # User accuracy
        self.assertTrue(
            np.allclose(sa.user_accuracy, np.array([0.97, 0.93, 0.97]))
        )
        # Producer accuracy
        self.assertTrue(
            np.allclose(
                sa.producer_accuracy,
                np.array([0.48063082, 0.99418868, 0.8969259]),
            )
        )
        # Overall accuracy
        self.assertEqual(sa.overall_accuracy, 0.9444167819481701)

    def test_map14(self):
        sa = MapSamples(class_area=[18000, 13500, 288000, 580500])
        sa.sample_size([0.6, 0.7, 0.9, 0.95], standard_error=0.01)

        # Sample size
        self.assertEqual(sa.n, 638)
        # Proportional allocation
        self.assertTrue(
            np.allclose(sa.prop_alloc, np.array([12, 9, 205, 413]))
        )
        # Proportional re-allocation
        self.assertTrue(
            np.allclose(sa.prop_realloc, np.array([50, 50, 165, 373]))
        )

    def test_maps4(self):
        error_matrix = np.array(
            [[48, 7, 0, 0], [13, 216, 0, 1], [1, 0, 49, 0], [3, 5, 0, 42]],
            dtype='float64',
        )

        sa = MapSamples(
            class_area=[47996, 228551, 13795, 3561, 293, 87],
            error_matrix=error_matrix,
        )
        sa.sample_size([0.01, 0.01, 0.0, 0.8, 0.0, 0.0], standard_error=0.005)

        # Sample size
        self.assertEqual(sa.n, 411)
        # Proportional allocation
        self.assertTrue(
            np.allclose(sa.prop_alloc, np.array([63, 300, 18, 4, 0, 0]))
        )
        # Proportional re-allocation
        self.assertTrue(
            np.allclose(sa.prop_realloc, np.array([50, 261, 50, 50, 0, 0]))
        )
