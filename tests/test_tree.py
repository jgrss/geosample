import unittest

import geopandas as gpd
import numpy as np
from shapely.geometry import Point

from geosample import QuadTree
from geosample.tree import GRTSFrame


class TestTree(unittest.TestCase):
    def test_query(self):
        data = np.array(
            [
                [-90, 40],
                [-91, 41],
                [-88, 38],
            ]
        )
        search_points = np.array([[-87, 39], [-88, 42]])
        geometry = [Point(xy) for xy in data]
        df = gpd.GeoDataFrame(
            data=range(data.shape[0]), geometry=geometry, crs='epsg:4326'
        )
        query_df = df.grts.query_points(search_points, k=1)
        self.assertEqual(len(query_df.index), search_points.shape[0])
        self.assertTrue(
            np.allclose(
                query_df.point_distance.values,
                np.array([1.41421356, 2.82842712]),
            )
        )
        self.assertTrue(
            np.allclose(np.array(query_df.index), np.array([2, 0]))
        )
        self.assertTrue(query_df.iloc[0].geometry == Point(-88, 38))
        self.assertTrue(query_df.iloc[1].geometry == Point(-90, 40))
