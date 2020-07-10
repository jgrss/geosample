from abc import ABC, abstractmethod
from collections import defaultdict, namedtuple

import numpy as np
import pandas as pd
import geopandas as gpd
from sklearn.cluster import KMeans
from shapely.geometry import Polygon


BBox = namedtuple('BBox', 'left bottom right top')


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

    """
    A class to generate a QuadTree

    Args:
        dataframe (DataFrame)
        force_square (Optional[bool])
    """

    def __init__(self, dataframe, force_square=True):

        super(QuadTree, self).__init__(dataframe)

        bounds_ = self.dataframe.total_bounds.flatten().tolist()
        bounds_names = self.bounds_to_tuple(bounds_)

        # Initiate the tree as the total bounds
        self.tree_bounds = [bounds_]
        self.tree_ids = ['0']

        if force_square:

            # Update the grid to force quadrants of equal length
            if self.min_qside == 'y':
                bounds_ = (bounds_names.left, bounds_names.top-abs(self.qx_len), bounds_names.right, bounds_names.top)
            else:
                bounds_ = (bounds_names.left, bounds_names.bottom, bounds_names.left+abs(self.qy_len), bounds_names.top)

            self.tree_bounds = [bounds_]

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

    @staticmethod
    def bounds_to_tuple(bounds):

        return BBox(left=bounds[0],
                    bottom=bounds[1],
                    right=bounds[2],
                    top=bounds[3])

    @property
    def min_qside(self):
        """Get the minimum quadrant side (y or x)"""
        return 'y' if self.qy_len < self.qx_len else 'x'

    @property
    def qy_len(self):
        """Get the quadrant latitudinal length"""
        bbox = self.bounds_to_tuple(self.tree_bounds[0])
        return bbox.top - bbox.bottom

    @property
    def qx_len(self):
        """Get the quadrant longitudinal length"""
        bbox = self.bounds_to_tuple(self.tree_bounds[0])
        return bbox.right - bbox.left

    @property
    def qmin(self):
        """Get the minimum quadrant length"""
        return self.qy_len if self.min_qside == 'y' else self.qx_len

    @property
    def qmax(self):
        """Get the maximum quadrant length"""
        return self.qy_len if self.min_qside == 'x' else self.qx_len

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
        Get counts of sample occurrences in each quadrant
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
        Splits a tree into quadrants

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

                if self.qmax <= max_length:
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

    def weight_grids(self, n_clusters=10, num_results=2):

        """
        Weights grids for sampling

        Args:
            n_clusters (Optional[int]): The number of clusters.
            num_results (Optional[int]): The number of result near cluster centers.

        Returns:
            ``geopandas.DataFrame``
        """

        qt_frame = self.to_frame()

        # Get coordinates
        X = np.c_[qt_frame.centroid.x.values, qt_frame.centroid.y.values]

        # Fit a KMeans
        kmeans = KMeans(n_clusters=n_clusters).fit(X)

        # Get the n nearest grids to the cluster centers
        for cluster_index in range(0, kmeans.cluster_centers_.shape[0]):

            bounds = (kmeans.cluster_centers_[cluster_index, 0],
                      kmeans.cluster_centers_[cluster_index, 1],
                      kmeans.cluster_centers_[cluster_index, 0],
                      kmeans.cluster_centers_[cluster_index, 1])

            sindex = qt_frame.sindex
            near_clusters = list(sindex.nearest(bounds, num_results=num_results))

            # Duplicate the near grids
            qt_frame = pd.concat((qt_frame, qt_frame.iloc[near_clusters]), axis=0)

        return qt_frame

    def sample(self, n=None, weight_by_clusters=False, random_state=None, **kwargs):

        """
        Samples from the hierarchical grid address using the
        Generalized Random Tessellation Stratified (GRTS) method

        Args:
            n (int): The target sample size.
            weight_by_clusters (Optional[bool])
            random_state (Optional[int])
            kwargs (Optional[dict]): Keyword arguments for ``self.weight_grids``.

        Returns:
            ``geopandas.GeoDataFrame``
        """

        if not isinstance(n, int):
            n = 1

        if isinstance(random_state, int):
            np.random.seed(random_state)

        # Sort by base 4 id
        if weight_by_clusters:
            df = self.weight_grids(**kwargs)
        else:
            df = self.to_frame()

        df = df.sort_values(by='id')

        npool = df.shape[0]
        interval = int(np.ceil(npool / n))

        # Get a random starting index
        start = np.random.randint(0, high=interval, size=1, dtype=int)[0]

        # Get the sample indices
        sample_indices = np.arange(start, npool, interval)

        # Get the random grids
        df_sample = df.iloc[sample_indices]

        sample_indices = []

        # Iterate over the selected grids,
        # get intersecting samples, and
        # select 1 sample within each grid.
        for row in df_sample.itertuples():

            # The grid bounds
            bbox = row.geometry.bounds

            # Points that intersect the current grid
            point_int = list(self.sindex.intersection(bbox))

            # Get one random point within the grid
            sample_indices.append(np.random.choice(point_int, size=1, replace=False)[0])

        # Get the random points
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
