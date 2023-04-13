import numpy as np
import scipy.stats as st


class MapSamples(object):
    """A class for sampling for probability-based map change."""

    def __init__(self, class_area, error_matrix=None, conf=0.95):

        self.class_area = class_area
        self.error_matrix = error_matrix
        self.conf = conf

        self.n = None
        self.prop_alloc = None
        self.prop_realloc = None

        if not isinstance(self.class_area, np.ndarray):
            self.class_area = np.array(self.class_area, dtype='float64')

        if isinstance(error_matrix, np.ndarray):

            if self.error_matrix.dtype != 'float64':
                self.error_matrix = np.float64(self.error_matrix)

    def test_13(self):
        """Test example from 2013 paper.

        Example:
            >>> sa = MapSamples([22353, 1122543, 610228])
            >>> sa.test_13()
        """

        self.error_matrix = np.array(
            [[97, 0, 3], [3, 279, 18], [2, 1, 97]], dtype='float64'
        )

        print('Class weights:', self.w)
        print('Prop. class area:', self.p_j)
        print('Stratified unbiased area estimate:', self.a_j)
        print('Standard error:', self.s_pj)
        # print('Standard error of the error-adjusted area:', sa.s_aj)
        print(
            'Standard error of the error-adjusted area conf.:', self.s_aj_conf
        )
        print('User accuracy:', self.user_accuracy)
        print('Producer accuracy:', self.producer_accuracy)
        print('Overall accuracy:', self.overall_accuracy)

    def test_14(self):
        """Test example from 2014 paper.

        Example:
            >>> sa = MapSamples([18000, 13500, 288000, 580500])
            >>> sa.test_14()
        """
        self.sample_size([0.6, 0.7, 0.9, 0.95], standard_error=0.01)
        print('Sample size:', self.n)
        print('Prop. allocation:', self.prop_alloc)
        print('Prop. re-allocation:', self.prop_realloc)

    def test_s4(self):
        """Test example from online PDF.

        Reference:
            https://github.com/beeoda/tutorials/blob/master/4_Estimation/S4_Methods_Estimation.pdf

        Example:
            >>> error_matrix = np.array([[48, 7, 0, 0], [13, 216, 0, 1], [1, 0, 49, 0], [3, 5, 0, 42]], dtype='float64')
            >>> sa = MapSamples([47996, 228551, 13795, 3561, 293, 87])
            >>> sa = MapSamples([47996, 228551, 13795, 3561], error_matrix=error_matrix)
        """

        self.sample_size(
            [0.01, 0.01, 0.0, 0.8, 0.0, 0.0], standard_error=0.005
        )
        print('Sample size:', self.n)
        print('Prop. allocation:', self.prop_alloc)
        print('Prop. re-allocation:', self.prop_realloc)

    def sample_size(
        self, users_accuracy, standard_error=0.01, interval=50, min_samples=0
    ):
        """Calculates the sample size given a target standard error.

        Args:
            users_accuracy (list): A list of class user accuracies.
            standard_error (Optional[float]): The target standard error.
            interval (Optional[int]): The re-allocation interval (e.g., 50, 75, 100).
            min_samples (Optional[int]): The minimum number of samples required to exclude from re-allocation.

        References:

            Olofsson et al. (2014) Good practices for estimating area and assessing
                accuracy of land change. Remote Sensing of Environment 148, 442-57.

                Section 5.1.1, page 53, Equation 13

        Example:
            # Example 5
            >>> sa = MapSamples()
            >>>
            >>> class_area = [18000, 13500, 288000, 580500]
            >>>
            >>> # Estimates of user's accuracy for each class
            >>> users_accuracy = [0.7, 0.6, 0.9, 0.95]
            >>>
            >>> sa.sample_size(class_area, users_accuracy, standard_error=0.01)
            >>> sa.n --> 641
        """
        if not isinstance(users_accuracy, np.ndarray):
            users_accuracy = np.array(users_accuracy, dtype='float64')

        # Estimated user's accuracy variance
        # Eq. 6
        v_u = users_accuracy * (1.0 - users_accuracy)

        # Stratum standard deviation
        # Cochran (1977)
        # Eq. 5.55
        s = np.sqrt(v_u)

        self.n = int(round(((self.w * s).sum() / standard_error) ** 2))
        self.prop_alloc = np.int64(self.n * self.w)
        self.prop_realloc = self.prop_alloc.copy()

        # Attempt to readjust
        attempts = 0
        while (
            self.prop_realloc[np.nonzero(self.prop_realloc)[0]].min()
            < interval
        ):

            adj1 = np.where(
                self.prop_realloc <= min_samples,
                0,
                interval - self.prop_realloc,
            )
            adj2 = np.where(
                (adj1 > min_samples) & (adj1 < interval),
                self.prop_realloc + adj1,
                self.prop_realloc,
            )
            diff = adj2.sum() - self.prop_realloc.sum()
            dist = (adj1 < 0).sum()

            adj3 = np.where(
                adj1 < 0, np.int64(self.prop_realloc - diff / dist), adj2
            )
            self.prop_realloc = np.where(
                (self.prop_alloc >= interval) & (adj3 < interval),
                interval,
                adj3,
            )

            if attempts > 10:
                break

            attempts += 1

        self.n = self.prop_realloc.sum()

    @property
    def user_accuracy(self):
        """Get the user's accuracy."""
        return np.diagonal(
            self.error_matrix_prop
        ) / self.error_matrix_prop.sum(axis=0)

    @property
    def producer_accuracy(self):
        """Get the producer's accuracy."""
        return np.diagonal(
            self.error_matrix_prop
        ) / self.error_matrix_prop.sum(axis=1)

    @property
    def overall_accuracy(self):
        """Get the overall accuracy."""
        return np.diagonal(self.error_matrix_prop).sum()

    @property
    def error_matrix_prop(self):
        return np.array(
            [
                self.w * (self.error_matrix[:, j] / self.n_j)
                for j in range(0, self.q)
            ],
            dtype='float64',
        )

    @property
    def q(self):
        """Get the number of q categories."""
        return self.error_matrix.shape[0]

    @property
    def p_j(self):
        """Get an estimator of the proportional area for each class, Olofsson
        et al.

        (2013), Eq. 1
        """
        return np.array(
            [
                (self.w * (self.error_matrix[:, j].flatten() / self.n_j)).sum()
                for j in range(0, self.q)
            ],
            dtype='float64',
        )

    @property
    def a_j(self):
        """Get an unbiased estimator of the total area, Olofsson et al.

        (2013), Eq. 2
        """
        return self.a_tot * self.p_j

    @property
    def n_j(self):
        """Get the column totals."""
        return self.error_matrix.sum(axis=1)

    @property
    def n_i(self):
        """Get the row totals."""
        return self.error_matrix.sum(axis=0)

    @property
    def w(self):
        """Get the class area proportions."""
        return self.class_area / self.a_tot

    @property
    def a_tot(self):
        """Get the total area of the map."""
        return self.class_area.sum()

    @property
    def s_pj(self):
        """Get the standard error of the estimated area proportion, Olofsson et
        al.

        (2013), Eq. 3
        """
        return np.array(
            [
                (
                    self.w**2
                    * (
                        (
                            (self.error_matrix[:, j].flatten() / self.n_j)
                            * (
                                1.0
                                - self.error_matrix[:, j].flatten() / self.n_j
                            )
                        )
                        / (self.n_j - 1.0)
                    )
                ).sum()
                ** 0.5
                for j in range(0, self.q)
            ],
            dtype='float64',
        )

    @property
    def s_aj(self):
        """Get the standard error of the error-adjusted estimated area,
        Olofsson et al.

        (2013), Eq. 4
        """
        return self.a_tot * self.s_pj

    @property
    def s_aj_conf(self):
        """Get the standard error of the error-adjusted estimated area
        confidence interval, Olofsson et al.

        (2013), Eq. 5
        """
        return self.s_aj * st.norm.interval(self.conf, loc=0, scale=1)[1]
