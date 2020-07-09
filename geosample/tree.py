from abc import ABC, abstractmethod
from collections import defaultdict

import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon


class TreeMixin(ABC):

    def __init__(self, dataframe):
        self.dataframe = dataframe
        self.sindex = self.dataframe.sindex

    @property
    def crs(self):
        """Get the GeoDataFrame CRS"""
        return self.dataframe.crs

    def create_poly(self, bounds):

        """
        Creates a Polygon geometry

        Args:
            bounds (tuple): Left, bottom, right, top.

        Returns:
            ``shapely.geometry.Polygon``
        """

        left, bottom, right, top = bounds

        return Polygon([(left, bottom),
                        (left, top),
                        (right, top),
                        (right, bottom),
                        (left, bottom)])

    @abstractmethod
    def to_geom(self):
        pass

    @abstractmethod
    def to_frame(self):
        pass


class QuadTree(TreeMixin):

    def __init__(self, dataframe):

        super(QuadTree, self).__init__(dataframe)

        # Initiate the tree as the total bounds
        self.tree_bounds = [self.dataframe.total_bounds.flatten().tolist()]
        self.tree_ids = ['0']

    def __iter__(self):
        return self

    def __next__(self):
        self.split()

    @property
    def nquads(self):
        """Get the number of quadrants in the tree"""
        return len(self.tree)

    @property
    def tree(self):
        """Get the quadrant tree geometry"""
        return self.to_geom()

    @property
    def qlen(self):
        """Get the length of a quadrant"""
        return self.tree_bounds[0][2] - self.tree_bounds[0][0]

    def to_geom(self):
        """Converts quadrant bounds to geometry"""
        return [self.create_poly(bbox) for bbox in self.tree_bounds]

    def to_frame(self):
        """Converts tree quadrants to a DataFrame"""
        return gpd.GeoDataFrame(data=self.tree_ids,
                                geometry=self.to_geom(),
                                crs=self.crs,
                                columns=['id'])

    @property
    def counts(self):

        """
        Counts sample occurrences in each quadrant
        """

        counts = defaultdict(int)

        for i, geom in zip(self.tree_ids, self.tree):
            point_int = list(self.sindex.intersection(geom.bounds))
            if point_int:
                counts[i] += len(point_int)

        return dict(counts)

    def count(self, qid):

        """
        Counts sample occurrences in a quadrant

        Args:
            qid (str): The quadrant id.
        """

        bbox = self.to_frame().query(f"id == '{qid}'").geometry.bounds.values.flatten().tolist()

        # Get points that intersect the quadrant
        point_int = list(self.sindex.intersection(bbox))

        return len(point_int) if point_int else 0

    def split(self):

        """
        Splits a tree

        1 | 3
        --|--
        0 | 2
        """

        new_tree_bounds = []
        new_tree_ids = []

        self.contains_null = False

        for qi, quad in enumerate(self.tree):

            left, bottom, right, top = quad.bounds
            xcenter = left + (right - left) / 2.0
            ycenter = top - (top - bottom) / 2.0

            quad_id = self.tree_ids[qi]

            for id_, bbox in zip([1, 3, 0, 2],
                                 [(left, ycenter, xcenter, top),
                                  (xcenter, ycenter, right, top),
                                  (left, bottom, xcenter, ycenter),
                                  (xcenter, bottom, right, ycenter)]):

                if list(self.sindex.intersection(bbox)):
                    new_tree_bounds.append(bbox)
                    new_tree_ids.append(quad_id + str(id_))
                else:
                    self.contains_null = True

        self.tree_bounds = new_tree_bounds
        self.tree_ids = new_tree_ids

        return self

    def split_recursive(self, max_samples=100, max_length=None, first_null=False):

        """
        Splits quadrants recursively

        Args:
            max_samples (Optional[int]): The maximum number of samples.
            max_length (Optional[float]): The maximum length of a quadrant side. Overrides ``max_samples``.
            first_null (Optional[bool]): Whether to break on the first null quadrant. Overrides ``max_samples``.
        """

        old_count = 1e9

        while True:

            self.split()
            max_count = self.counts[max(self.counts, key=self.counts.get)]

            if isinstance(max_length, float) or isinstance(max_length, int):

                if self.qlen <= max_length:
                    break

            elif first_null:

                if self.contains_null:
                    break

            else:

                if max_count <= max_samples:
                    break

                if max_count == old_count:
                    break

            old_count = max_count

    def sample(self, n=None, random_state=None):

        """
        Samples from the hierarchical grid address using the
        Generalized Random Tessellation Stratified (GRTS) method

        Args:
            n (int): The target sample size.
            random_state (Optional[int])

        Returns:
            ``geopandas.GeoDataFrame``
        """

        if not isinstance(n, int):
            n = 1

        if isinstance(random_state, int):
            np.random.seed(random_state)

        df = self.to_frame().sort_values(by='id')

        npool = df.shape[0]
        interval = int(np.ceil(npool / n))

        start = np.random.randint(0, high=interval, size=1, dtype=int)[0]

        sample_indices = np.arange(start, npool, interval)

        df_sample = df.iloc[sample_indices]

        sample_indices = []

        for row in df_sample.itertuples():
            bbox = row.geometry.bounds

            point_int = list(self.sindex.intersection(bbox))

            sample_indices.append(np.random.choice(point_int, size=1, replace=False)[0])

        return self.dataframe.iloc[sample_indices]


class Rtree(TreeMixin):

    def __init__(self, dataframe):
        super(Rtree, self).__init__(dataframe)

    def __len__(self):
        for group_idx, indices, bbox in self.sindex.leaves():
            n = len(indices)
            break

        return n

    @property
    def nleaves(self):
        return len(self.sindex.leaves())

    def to_geom(self):
        """Converts leaves to geometry"""
        return [self.create_poly(bbox) for group_idx, indices, bbox in self.sindex.leaves()]

    def to_frame(self):
        """Converts leaves to a DataFrame"""
        return gpd.GeoDataFrame(data=range(0, self.nleaves),
                                geometry=self.to_geom(),
                                crs=self.crs,
                                columns=['id'])
