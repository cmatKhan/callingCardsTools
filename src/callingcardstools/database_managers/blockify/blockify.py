# import ptpdb
from functools import reduce
import numpy as np
import numpy.ma as ma
import pandas as pd
# from profilehooks import profile, timecall
from scipy.stats import poisson
import sys
import warnings


__all__ = ['blockify','Algorithm', 'BayesianBlocks', 'OptimalPartitioning', 'PELT', 'FPOP']

def blockify(t, x=None, sigma=None, algorithm='bayesian_blocks', **kwargs):
    """_summary_

    Args:
        t (_type_): _description_
        x (_type_, optional): _description_. Defaults to None.
        sigma (_type_, optional): _description_. Defaults to None.
        algorithm (str, optional): _description_. Defaults to 'bayesian_blocks'.

    Raises:
        ValueError: _description_

    Returns:
        _type_: _description_
    """
    # Compute optimal segmentation of data with Scargle's Bayesian Blocks

    # This is a flexible implementation of the Bayesian Blocks algorithm
    # described in Scargle 2012 [1]_.

    # Parameters
    # ----------
    # t : array_like
    #     data times (one dimensional, length N)
    # x : array_like (optional)
    #     data values
    # sigma : array_like or float (optional)
    #     data errors
    # algorithm : str or object
    #     the algorithm function to use for the model.
    #     If a string, the following options are supported:

    #     - 'events' : binned or unbinned event data.  Arguments are ``gamma``,
    #       which gives the slope of the prior on the number of bins, or
    #       ``ncp_prior``, which is :math:`-\ln({\tt gamma})`.
    #     - 'regular_events' : non-overlapping events measured at multiples of a
    #       fundamental tick rate, ``dt``, which must be specified as an
    #       additional argument.  Extra arguments are ``p0``, which gives the
    #       false alarm probability to compute the prior, or ``gamma``, which
    #       gives the slope of the prior on the number of bins, or ``ncp_prior``,
    #       which is :math:`-\ln({\tt gamma})`.
    #     - 'measures' : algorithm for a measured sequence with Gaussian errors.
    #       Extra arguments are ``p0``, which gives the false alarm probability
    #       to compute the prior, or ``gamma``, which gives the slope of the
    #       prior on the number of bins, or ``ncp_prior``, which is
    #       :math:`-\ln({\tt gamma})`.

    #     In all three cases, if more than one of ``p0``, ``gamma``, and
    #     ``ncp_prior`` is chosen, ``ncp_prior`` takes precedence over ``gamma``
    #     which takes precedence over ``p0``.

    #     Alternatively, the algorithm parameter can be an instance of
    #     :class:`Algorithm` or a subclass thereof.

    # **kwargs :
    #     any additional keyword arguments will be passed to the specified
    #     :class:`Algorithm` derived class.

    # Returns
    # -------
    # edges : ndarray
    #     array containing the (N+1) edges defining the N bins

    # Examples
    # --------
    # Event data:

    # >>> t = np.random.normal(size=100)
    # >>> edges = bayesian_blocks(t, algorithm='events', p0=0.01)

    # Event data with repeats:

    # >>> t = np.random.normal(size=100)
    # >>> t[80:] = t[:20]
    # >>> edges = bayesian_blocks(t, algorithm='events', p0=0.01)

    # Regular event data:

    # >>> dt = 0.05
    # >>> t = dt * np.arange(1000)
    # >>> x = np.zeros(len(t))
    # >>> x[np.random.randint(0, len(t), len(t) // 10)] = 1
    # >>> edges = bayesian_blocks(t, x, algorithm='regular_events', dt=dt)

    # Measured point data with errors:

    # >>> t = 100 * np.random.random(100)
    # >>> x = np.exp(-0.5 * (t - 50) ** 2)
    # >>> sigma = 0.1
    # >>> x_obs = np.random.normal(x, sigma)
    # >>> edges = bayesian_blocks(t, x_obs, sigma, algorithm='measures')

    # References
    # ----------
    # .. [1] Scargle, J et al. (2012)
    #    http://adsabs.harvard.edu/abs/2012arXiv1207.5578S

    # See Also
    # --------
    # astropy.stats.histogram : compute a histogram using bayesian blocks
    
    ALGORITHM_DICT = {'bayesian_blocks': BayesianBlocks,
                      'OP': OptimalPartitioning,
                      'PELT': PELT,
                      'FPOP': FPOP}
    algorithm = ALGORITHM_DICT.get(algorithm, algorithm)

    if type(algorithm) is type and issubclass(algorithm, Algorithm):
        alg = algorithm(**kwargs)
    elif isinstance(algorithm, Algorithm):
        alg = algorithm
    else:
        raise ValueError("algorithm parameter not understood")

    return alg.segment(t, x, sigma)


class Algorithm(object):
    """Base class for bayesian blocks algorithm functions

    Derived classes should overload the following method:

    ``algorithm(self, **kwargs)``:
      Compute the algorithm given a set of named arguments.
      Arguments accepted by algorithm must be among ``[T_k, N_k, a_k, b_k, c_k]``
      (See [1]_ for details on the meaning of these parameters).

    Additionally, other methods may be overloaded as well:

    ``__init__(self, **kwargs)``:
      Initialize the algorithm function with any parameters beyond the normal
      ``p0`` and ``gamma``.

    ``validate_input(self, t, x, sigma)``:
      Enable specific checks of the input data (``t``, ``x``, ``sigma``)
      to be performed prior to the fit.

    ``compute_ncp_prior(self, N)``: If ``ncp_prior`` is not defined explicitly,
      this function is called in order to define it before fitting. This may be
      calculated from ``gamma``, ``p0``, or whatever method you choose.

    ``p0_prior(self, N)``:
      Specify the form of the prior given the false-alarm probability ``p0``
      (See [1]_ for details).

    For examples of implemented algorithm functions, see :class:`Events`,
    :class:`RegularEvents`, and :class:`PointMeasures`.

    References
    ----------
    .. [1] Scargle, J et al. (2012)
       http://adsabs.harvard.edu/abs/2012arXiv1207.5578S
    """
    def __init__(self, p0=0.05, gamma=None, ncp_prior=None):
        self.p0 = p0
        self.gamma = gamma
        self.ncp_prior = ncp_prior

    def validate_input(self, t, x=None, sigma=None):
        """Validate inputs to the model.

        Parameters
        ----------
        t : array_like
            times of observations
        x : array_like (optional)
            values observed at each time
        sigma : float or array_like (optional)
            errors in values x

        Returns
        -------
        t, x, sigma : array_like, float or None
            validated and perhaps modified versions of inputs
        """
        # validate array input
        t = np.asarray(t, dtype=float)
        if x is not None:
            x = np.asarray(x)
        if sigma is not None:
            sigma = np.asarray(sigma)

        # find unique values of t
        t = np.array(t)
        if t.ndim != 1:
            raise ValueError("t must be a one-dimensional array")
        unq_t, unq_ind, unq_inv = np.unique(t, return_index=True,
                                            return_inverse=True)

        # if x is not specified, x will be counts at each time
        if x is None:
            if sigma is not None:
                raise ValueError("If sigma is specified, x must be specified")
            else:
                sigma = 1

            if len(unq_t) == len(t):
                x = np.ones_like(t)
            else:
                x = np.bincount(unq_inv)

            t = unq_t

        # if x is specified, then we need to simultaneously sort t and x
        else:
            # TODO: allow broadcasted x?
            x = np.asarray(x)
            if x.shape not in [(), (1,), (t.size,)]:
                raise ValueError("x does not match shape of t")
            x += np.zeros_like(t)

            if len(unq_t) != len(t):
                raise ValueError("Repeated values in t not supported when "
                                 "x is specified")
            t = unq_t
            x = x[unq_ind]

        # verify the given sigma value
        if sigma is None:
            sigma = 1
        else:
            sigma = np.asarray(sigma)
            if sigma.shape not in [(), (1,), (t.size,)]:
                raise ValueError('sigma does not match the shape of x')

        return t, x, sigma

    def fitness(self, **kwargs):
        raise NotImplementedError()

    def p0_prior(self, N):
        """
        Empirical prior, parametrized by the false alarm probability ``p0``
        See  eq. 21 in Scargle (2012)

        Note that there was an error in this equation in the original Scargle
        paper (the "log" was missing). The following corrected form is taken
        from http://arxiv.org/abs/1304.2818
        """
        return 4 - np.log(73.53 * self.p0 * (N ** -0.478))

    # the algorithm_args property will return the list of arguments accepted by
    # the method algorithm().  This allows more efficient computation below.
    @property
    def _algorithm_args(self):
        return signature(self.algorithm).parameters.keys()

    def compute_ncp_prior(self, N):
        """
        If ``ncp_prior`` is not explicitly defined, compute it from ``gamma``
        or ``p0``.
        """
        if self.ncp_prior is not None:
            return self.ncp_prior
        elif self.gamma is not None:
            return -np.log(self.gamma)
        elif self.p0 is not None:
            return self.p0_prior(N)
        else:
            raise ValueError("``ncp_prior`` is not defined, and cannot compute "
                             "it as neither ``gamma`` nor ``p0`` is defined.")
    # @timecall
    def segment(self, t, x=None, sigma=None):
        raise NotImplementedError()


class OptimalPartitioning(Algorithm):
    """Bayesian blocks algorithm for regular events

    This is for data which has a fundamental "tick" length, so that all
    measured values are multiples of this tick length.  In each tick, there
    are either zero or one counts.

    Parameters
    ----------
    dt : float
        tick rate for data
    p0 : float (optional)
        False alarm probability, used to compute the prior on :math:`N_{\rm
        blocks}` (see eq. 21 of Scargle 2012). If gamma is specified, p0 is
        ignored.
    ncp_prior : float (optional)
        If specified, use the value of ``ncp_prior`` to compute the prior as
        above, using the definition :math:`{\tt ncp\_prior} = -\ln({\tt
        gamma})`.  If ``ncp_prior`` is specified, ``gamma`` and ``p0`` are
        ignored.
    """
    def __init__(self, p0=0.05, gamma=None, ncp_prior=None):
        if p0 is not None and gamma is None and ncp_prior is None:
            warnings.warn('p0 does not seem to accurately represent the false '
                          'positive rate for event data. It is highly '
                          'recommended that you run random trials on signal-'
                          'free noise to calibrate ncp_prior to achieve a '
                          'desired false positive rate.')
        super(OptimalPartitioning, self).__init__(p0, gamma, ncp_prior)

    def fitness(self, T_k, N_k):
        # Negative log of the Poisson maximum likelihood, given T_k and N_k
        return N_k * (np.log(N_k) - np.log(T_k))

    # @timecall
    def segment(self, t, x=None, sigma=None):
        """Fit the Bayesian Blocks model given the specified algorithm function.

        Parameters
        ----------
        t : array_like
            data times (one dimensional, length N)
        x : array_like (optional)
            data values
        sigma : array_like or float (optional)
            data errors

        Returns
        -------
        edges : ndarray
            array containing the (M+1) edges defining the M optimal bins
        """

        t, x, sigma = self.validate_input(t, x, sigma)

        # # compute values needed for computation, below
        # if 'a_k' in self._algorithm_args:
        #     ak_raw = np.ones_like(x) / sigma ** 2
        # if 'b_k' in self._algorithm_args:
        #     bk_raw = x / sigma ** 2
        # if 'c_k' in self._algorithm_args:
        #     ck_raw = x * x / sigma ** 2

        # create length-(N + 1) array of cell edges
        edges = np.concatenate([t[:1],
                                0.5 * (t[1:] + t[:-1]),
                                t[-1:]])
        block_length = t[-1] - edges

        # arrays to store the best configuration
        N = len(t)
        print(N)
        best = np.zeros(N, dtype=float)
        last = np.zeros(N, dtype=int)

        # Compute ncp_prior if not defined
        if self.ncp_prior is None:
            ncp_prior = self.compute_ncp_prior(N)
        # ----------------------------------------------------------------
        # Start with first data cell; add one cell at each iteration
        # ----------------------------------------------------------------
        kwds = {}
        for R in range(N):
            # Compute fit_vec : algorithm of putative last block (end at R)

            # T_k: width/duration of each block
            kwds['T_k'] = block_length[:R + 1] - block_length[R + 1]

            # N_k: number of elements in each block
            kwds['N_k'] = np.cumsum(x[:R + 1][::-1])[::-1]

            # evaluate algorithm function
            fit_vec = self.fitness(**kwds)
            # print(np.exp(-fit_vec))
            A_R = fit_vec - ncp_prior
            # print(A_R)
            A_R[1:] += best[:R]

            i_max = np.argmax(A_R)
            last[R] = i_max
            best[R] = A_R[i_max]
            # print("best:")
            # print(best)
            # print("fit_vec:")
            # print(fit_vec)
            # print("last:")
            # print(last)
            # print()
        print(best[R])
        # ----------------------------------------------------------------
        # Now find changepoints by iteratively peeling off the last block
        # ----------------------------------------------------------------
        change_points = np.zeros(N, dtype=int)
        i_cp = N
        ind = N
        while True:
            i_cp -= 1
            change_points[i_cp] = ind
            if ind == 0:
                break
            ind = last[ind - 1]
        change_points = change_points[i_cp:]
        return edges[change_points]


class BayesianBlocks(OptimalPartitioning):
    r"""Bayesian blocks algorithm for binned or unbinned events

    Parameters
    ----------
    p0 : float (optional)
        False alarm probability, used to compute the prior on
        :math:`N_{\rm blocks}` (see eq. 21 of Scargle 2012). For the Events
        type data, ``p0`` does not seem to be an accurate representation of the
        actual false alarm probability. If you are using this algorithm function
        for a triggering type condition, it is recommended that you run
        statistical trials on signal-free noise to determine an appropriate
        value of ``gamma`` or ``ncp_prior`` to use for a desired false alarm
        rate.
    gamma : float (optional)
        If specified, then use this gamma to compute the general prior form,
        :math:`p \sim {\tt gamma}^{N_{\rm blocks}}`.  If gamma is specified, p0
        is ignored.
    ncp_prior : float (optional)
        If specified, use the value of ``ncp_prior`` to compute the prior as
        above, using the definition :math:`{\tt ncp\_prior} = -\ln({\tt
        gamma})`.
        If ``ncp_prior`` is specified, ``gamma`` and ``p0`` is ignored.
    """
    def __init__(self, p0=0.05, gamma=None, ncp_prior=None):
        super(BayesianBlocks, self).__init__(p0, gamma, ncp_prior)

    def fitness(self, N_k, T_k):
        # eq. 19 from Scargle 2012
        return N_k * (np.log(N_k) - np.log(T_k))
        # return poisson.logpmf(N_k, N_k/T_k)

    # def validate_input(self, t, x, sigma):
    #     t, x, sigma = super(BayesianBlocks, self).validate_input(t, x, sigma)
    #     if x is not None and np.any(x % 1 > 0):
    #         raise ValueError("x must be integer counts for algorithm='events'")
    #     return t, x, sigma


class PELT(Algorithm):
    r"""Bayesian blocks algorithm for point measures

    Parameters
    ----------
    p0 : float (optional)
        False alarm probability, used to compute the prior on :math:`N_{\rm
        blocks}` (see eq. 21 of Scargle 2012). If gamma is specified, p0 is
        ignored.
    ncp_prior : float (optional)
        If specified, use the value of ``ncp_prior`` to compute the prior as
        above, using the definition :math:`{\tt ncp\_prior} = -\ln({\tt
        gamma})`.  If ``ncp_prior`` is specified, ``gamma`` and ``p0`` are
        ignored.
    """
    def __init__(self, p0=0.05, gamma=None, ncp_prior=None):
        super(PELT, self).__init__(p0, gamma, ncp_prior)

    def fitness(self, N_k, T_k):
        # eq. 19 from Scargle 2012
        return -N_k * (np.log(N_k) - np.log(T_k))
        # return -poisson.logpmf(N_k, N_k / T_k)

    # @timecall
    def segment(self, t, x=None, sigma=None):
        """Fit the Bayesian Blocks model given the specified algorithm function.

        Parameters
        ----------
        t : array_like
            data times (one dimensional, length N)
        x : array_like (optional)
            data values
        sigma : array_like or float (optional)
            data errors

        Returns
        -------
        edges : ndarray
            array containing the (M+1) edges defining the M optimal bins
        """

        t, x, sigma = self.validate_input(t, x, sigma)
        # # compute values needed for computation, below
        # if 'a_k' in self._algorithm_args:
        #     ak_raw = np.ones_like(x) / sigma ** 2
        # if 'b_k' in self._algorithm_args:
        #     bk_raw = x / sigma ** 2
        # if 'c_k' in self._algorithm_args:
        #     ck_raw = x * x / sigma ** 2

        # create length-(N + 1) array of cell edges
        edges = np.concatenate([t[:1],
                                0.5 * (t[1:] + t[:-1]),
                                t[-1:]])
        block_length = t[-1] - edges

        # arrays to store the best configuration
        N = len(t)
        # print(N)
        best = np.zeros(N + 1, dtype=float)
        last = np.zeros(N, dtype=int)

        # Compute ncp_prior if not defined
        if self.ncp_prior is None:
            ncp_prior = self.compute_ncp_prior(N)
        K = 0
        unpruned = np.array([], dtype = np.int)
        lastPruned  = -1
        # ----------------------------------------------------------------
        # Start with first data cell; add one cell at each iteration
        # ----------------------------------------------------------------
        for R in range(N):
            # Consider everything unpruned until proven otherwise
            unpruned = np.concatenate([unpruned, np.array([R])])
            # T_k: width/duration of each block
            T_k = block_length[unpruned] - block_length[R + 1]
            # N_k: number of elements in each block
            N_k = np.cumsum(x[unpruned][::-1])[::-1]
            # Compute fit_vec : algorithm of putative last block (end at R)
            fit_vec = self.fitness(N_k, T_k)
            A_R = fit_vec + ncp_prior# + best[lastPruned]
            # print("best: ", best[:R])
            # A_R[1:] += best[lastPruned + 1:R] # before
            A_R += best[lastPruned + 1:R + 1] # after
            i_min = np.argmin(A_R)#optint[unpruned])
            last[R] = unpruned[i_min]
            # best[R] = A_R[i_min] # before
            best[R + 1] = A_R[i_min] # after
            # peltFitness = best[lastPruned + 1:R] + fit_vec[1:] # before
            peltFitness = best[lastPruned + 1:R + 1] + fit_vec # after
            # Pruning step; only applies if N_k > 2
            candidates = (N_k[:-1] >= 2).nonzero()[0]
            if candidates.size:
                pruned = (peltFitness + K > best[R+1]).nonzero()[0]
                if pruned.size:
                    case = 3
                    # Case 1: take the largest index and prune everything before that
                    i_max = np.max(pruned)
                    if case == 1:
                        pass
                    # Case 2: find the largest continuous value in pruned
                    elif case == 2:
                        i_max = 0
                        if len(pruned) > 1:
                            while i_max != len(pruned) - 1:
                                if pruned[i_max + 1] - pruned[i_max] == 1:
                                    i_max += 1
                                else:
                                    break
                    # Case 3: wait until all indexes from 0 to i_max satisfy the pruning condition; then prune i_max
                    elif case == 3:
                        if np.array_equal(pruned, np.arange(i_max + 1)):
                            i_max = np.max(pruned)
                        else:
                            i_max = -1
                    if i_max >= 0:
                        t_max = unpruned[i_max]
                        unpruned = np.delete(unpruned, np.arange(i_max + 1))
                        lastPruned = t_max
        # print(best[R+1])
        # ----------------------------------------------------------------
        # Now find changepoints by iteratively peeling off the last block
        # ----------------------------------------------------------------
        change_points = np.zeros(N, dtype=int)
        i_cp = N
        ind = N
        while True:
            i_cp -= 1
            change_points[i_cp] = ind
            if ind == 0:
                break
            ind = last[ind - 1]
        change_points = change_points[i_cp:]
        return edges[change_points]


class FPOP(Algorithm):
    """Bayesian blocks algorithm for point measures

    Parameters
    ----------
    p0 : float (optional)
        False alarm probability, used to compute the prior on :math:`N_{\rm
        blocks}` (see eq. 21 of Scargle 2012). If gamma is specified, p0 is
        ignored.
    ncp_prior : float (optional)
        If specified, use the value of ``ncp_prior`` to compute the prior as
        above, using the definition :math:`{\tt ncp\_prior} = -\ln({\tt
        gamma})`.  If ``ncp_prior`` is specified, ``gamma`` and ``p0`` are
        ignored.
    """
    def __init__(self, p0=0.05, gamma=None, ncp_prior=None):
        super(FPOP, self).__init__(p0, gamma, ncp_prior)

    def fitness(self, N_k, T_k):
        # eq. 19 from Scargle 2012
        return -N_k * (np.log(N_k) - np.log(T_k))

    def validate_input(self, t, x, sigma):
        if x is None:
            raise ValueError("x must be specified for point measures")
        return super(FPOP, self).validate_input(t, x, sigma)


# this was in the original file -- not sure exactly why...it seems like this 
# gets called from other scripts
#
# if __name__ == "__main__":
#     expFile = pd.read_table(sys.argv[1], header = None)
#     chroms = reduce(lambda l, x: l if x in l else l+[x], expFile[0], [])
#  
#     for chrom in chroms:
#         blocks = blockify(expFile[expFile[0] == chrom][1], algorithm="bayesian_blocks")
#         # print(blocks)