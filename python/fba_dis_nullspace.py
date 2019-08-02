import argparse
import sys
from multiprocessing import Pool
from collections import OrderedDict
import pandas as pd
from sklearn.decomposition import PCA
from scipy.linalg import null_space, lstsq
import numpy as np
from psamm.lpsolver.cplex import Solver
from psamm.fluxanalysis import FluxBalanceProblem, FluxBalanceError
from psamm.datasource.native import ModelReader
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa


class stepError(Exception):
    def __init__(self, *args, **kwargs):
        super(stepError, self).__init__(*args, **kwargs)


class sampler(object):
    def __init__(self, metabolic_model, solver, fix=None,
                 objective=None, objective_scale=1, nproc=1, epsilon=1e-5):
        self._model = metabolic_model
        self._solver = solver
        # the original stoichiometric matrix
        self._sm = np.array(self._model.matrix)
        # the null space of stoichiometric matrix
        self._ns = null_space(self._sm)
        self._compounds = sorted(self._model.compounds)
        self._reactions = sorted(self._model.reactions)
        self._epsilon = epsilon
        self._p = FluxBalanceProblem(self._model, self._solver)
        self._objective = objective
        self._objective_scale = objective_scale
        self._nproc = nproc
        self._fix = dict()
        # set value for fixed reactions
        if fix is not None:
            fix = fix.split(';')
            for i in fix:
                rxn, value = i.split(',')
                self._fix[rxn] = float(value)
                self._p.prob.add_linear_constraints(
                    self._p.get_flux_var(rxn) == float(value)
                )
        # set objective reactions
        if self._objective is not None:
            if not self._model.has_reaction(self._objective):
                raise RuntimeError('No reaction %s in model!' % objective)
            # get the maximium objective value of this model
            self._p.maximize(self._objective)
            self._objective_value = self._p.get_flux(self._objective)
            # force the model to product max objective flux
            self._p.prob.add_linear_constraints(self._p.get_flux_var(
                self._objective) >= self._objective_value * objective_scale)
        self._set_limits()

    def _set_limits(self):
        """Set the upper and lower bound of reaction flux"""
        self._upper = OrderedDict()
        self._lower = OrderedDict()
        for rxn in self._reactions:
            if rxn in self._fix:
                self._lower[rxn] = self._fix[rxn]
                self._upper[rxn] = self._fix[rxn]
            else:
                l, u = self._model.limits[rxn]
                self._lower[rxn] = float(l)
                self._upper[rxn] = float(u)
        # reset limits of objective reaction
        if self._objective is not None and self._objective not in self._fix:
            self._upper[self._objective] = self._objective_value
            self._lower[self._objective] = (self._objective_value
                                            * self._objective_scale)

    def _get_projection(self, x, check=True):
        """obtain the FBA result in null space"""
        projection = lstsq(self._ns, x)
        if projection[1] > self._epsilon and check:
            raise RuntimeError('Failed to project point into null space!')
        return projection[0]

    def _get_fluxes(self):
        """obtain the FBA result"""
        x = list()
        for rxn in self._reactions:
            x.append(self._p.get_flux(rxn))
        x = np.array(x)
        return x

    def set_warmup(self):
        """Set up warmup points for Monte Carlo sampling."""
        self._warmup = list()
        self._warmup_flux = list()
        print('Setting up warmup points...')
        for i in self._reactions:
            if self._upper[i] - self._lower[i] < self._epsilon:
                # skip fixed reaction
                continue
            try:  # maximize the flux of reaction i
                self._p.maximize(i)
                # store warmup points based on non-zero reacions only
                if np.abs(self._p.get_flux(i)) > self._epsilon:
                    fluxes = self._get_fluxes()
                    self._warmup_flux.append(fluxes)
                    self._warmup.append(self._get_projection(fluxes))
            except FluxBalanceError:
                pass
            if self._model.is_reversible(i):
                # maximize flux of reaction i in reverse direction
                try:
                    self._p.maximize({i: -1})
                    # store warmup points based on effective reacions only
                    if np.abs(self._p.get_flux(i)) > self._epsilon:
                        fluxes = self._get_fluxes()
                        self._warmup_flux.append(fluxes)
                        self._warmup.append(self._get_projection(fluxes))
                except FluxBalanceError:
                    pass
        self._warmup = np.array(self._warmup)
        if len(self._warmup) <= 1:
            raise RuntimeError("Can't get solutions based on current "
                               "flux limitations! Please check your model.")
        # maintain unrelated warmup points only
        self._warmup = non_redundant(self._warmup)
        self._warmup_flux = np.array(self._warmup_flux)

    @property
    def warmup_flux(self):
        """The reaction fluxes of warmup points"""
        return self._warmup_flux

    @property
    def p(self):
        """The underlying :class:`FluxBalanceProblem`"""
        return self._p


class optGp(sampler):
    def __init__(self, *args, **kwargs):
        super(optGp, self).__init__(*args, **kwargs)

    def sample(self, nsample, k=100, maxtry=1):
        """Artificial Centering Hit-and-Run functioin"""
        print('Doing artificial centering hit-and-run')
        # number of warmup points
        nwarm = len(self._warmup)
        maxiter = nsample * k // self._nproc + 1
        upper = np.array([i for i in self._upper.values()])
        lower = np.array([i for i in self._lower.values()])
        if self._nproc == 1:
            m = np.random.choice(nwarm)
            x = one_chain(m, self._warmup, maxiter, self._ns,
                          upper, lower, self._epsilon,
                          k, maxtry)
            return pd.DataFrame(x, columns=self._reactions)

        elif self._nproc > 1:
            pool = Pool(self._nproc)
            tasks = ((one_chain,
                      (m, self._warmup, maxiter, self._ns,
                       upper, lower, self._epsilon, k, maxtry))
                     for m in np.random.choice(nwarm, self._nproc))
            x = pool.map(parallel_worker, tasks)
            pool.close()
            pool.join()
            return pd.DataFrame(np.vstack(x),
                                columns=self._reactions).round(
                                    int(-np.log10(self._epsilon)))


def non_redundant(x, threshold=0.95):
    """Run pairwise Pearson correlation to find non-redundant rows"""
    corr = np.corrcoef(x)
    corr = np.tril(corr, -1)
    index = (corr > threshold).any(axis=1)
    return x[index]


def parallel_worker(tasks):
    """Send eclapsed parameters to function"""
    func, params = tasks
    return func(*params)


def one_chain(m, warmup, maxiter, ns, upper, lower, epsilon, k, maxtry):
    """Run one ACHR chain"""
    print('Start on point %i' % m)
    x = list()
    nwarm = len(warmup)
    npoints = nwarm
    # get the center point
    s = warmup.mean(axis=0)
    # set up starting point
    # pull back a bit to avoid stuck
    xm = (warmup[m] - s) * 0.9 + s
    xm_flux = ns.dot(xm)
    for niter in range(maxiter):
        success = False
        # wait until a successful move
        while True:
            for ntry in range(maxtry):
                n = np.random.choice(nwarm)
                xn = warmup[n]
                direction = xn - s
                direction_flux = ns.dot(direction)
                try:
                    xm, xm_flux = one_step(xm, xm_flux, direction,
                                           direction_flux, ns,
                                           upper, lower, epsilon)
                    success = True  # successfully got the next point
                    break
                except stepError:
                    pass
            if success:
                break
            # too many failures, move xm to new position
            # set up starting point
            # pull back to center to avoid stuck
            xm = s
            xm_flux = ns.dot(xm)
        npoints += 1
        # recalculate the center
        s = (s * npoints + xm) / (npoints + 1)
        # output only one point every k iter
        if (niter + 1) % k == 0:
            # add new point
            x.append(xm_flux)
            if (niter + 1) // k % 100 == 0:
                print((niter + 1) // k)
                sys.stdout.flush()
    print('Point %i is done...' % m)
    return np.array(x)


def one_step(xm, xm_flux, direction, direction_flux,
             ns, upper, lower, epsilon):
    """Move one step further on one chain"""
    index = np.abs(direction_flux) > epsilon
    up_zero = np.abs(upper - xm_flux)[index] < epsilon
    down_zero = np.abs(lower - xm_flux)[index] < epsilon
    scale_upper = (upper
                   - xm_flux)[index] / direction_flux[index]
    scale_upper[up_zero] = 0
    scale_lower = (lower
                   - xm_flux)[index] / direction_flux[index]
    scale_lower[down_zero] = 0
    scale = np.array([scale_upper, scale_lower])
    up = scale.max(axis=0).min()
    down = scale.min(axis=0).max()
    if up - down < epsilon:
        # can't move
        raise stepError('Cannot move!')
    # get the new point
    step = np.random.uniform(down, up)
    new = xm + step * direction
    new_flux = ns.dot(new)
    while (np.sum(new_flux - upper > epsilon) > 0
            or np.sum(new_flux - lower < -epsilon) > 0):
        # new point is out of boundary
        step = np.random.sample() * step
        new = xm + step * direction
        new_flux = ns.dot(new)
    # got acceptable result
    return new, new_flux


class ACHR(sampler):
    def __init__(self, *args, **kwargs):
        super(ACHR, self).__init__(*args, **kwargs)

    def _random_optimize(self):
        # new optimize project
        optimize = dict()
        for i in self._reactions:
            limit_range = self._upper[i] - self._lower[i]
            # skip fixed reaction
            if limit_range < self._epsilon:
                continue
            # randomly set optimize weight
            if self._upper[i] * self._lower[i] < 0:
                optimize[i] = 2 * np.random.sample() - 1.0
            else:
                optimize[i] = np.random.sample()
        self._p.maximize(optimize)
        xm_flux = self._get_fluxes()
        xm = self._get_projection(xm_flux)
        return xm, xm_flux

    def sample(self, nsample):
        """Artificial Centering Hit-and-Run functioin"""
        print('Doing artificial centering hit-and-run')
        points = [row for row in self._warmup]
        points_flux = [row for row in self._warmup_flux]
        npoint = len(points)
        nwarm = len(self._warmup)
        # center of warmup points
        s = self._warmup.mean(axis=0)
        # starting point
        # xm = points[npoint - 1]
        xm = s
        xm_flux = points_flux[npoint - 1]
        upper = np.array([i for i in self._upper.values()])
        lower = np.array([i for i in self._lower.values()])
        while npoint - nwarm < nsample:
            n = np.random.choice(npoint)
            xn = points[n]
            direction = xn - s
            direction_flux = self._ns.dot(direction)
            try:
                xm, xm_flux = one_step(xm, xm_flux, direction, direction_flux,
                                       self._ns, upper, lower, self._epsilon)
                points.append(xm)
                points_flux.append(xm_flux)
                # recalculate center
                s = (s * npoint + xm) / (npoint + 1)
                npoint += 1
                if (npoint - nwarm) % 100 == 0:
                    print('%i/%i' % (npoint - nwarm, nsample))
                    sys.stdout.flush()
            except stepError:
                # print('stuck')
                xm = s
                xm_flux = self._ns.dot(xm)
                pass
        return pd.DataFrame(points_flux[nwarm:],
                            columns=self._reactions).round(
                                int(-np.log10(self._epsilon)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            'Get the distribution of reaction fluxes through '
            'Artificial Centering Hit-and-Run algorithm.'
        )
    )
    parser.add_argument('--model', help='path to metabolic model',
                        required=True)
    parser.add_argument('-s', '--samples', type=int, default=10000,
                        help='number of samples to get (default: 10000)')
    parser.add_argument('--objective',
                        help=('reaction to use as objective, '
                              'type "biomass" for biomass reaction'))
    parser.add_argument('--threshold', type=float, default=1.0,
                        help=('ratio of maximum objective flux to maintain '
                              '(0.0~1.0, default: 1.0)'))
    parser.add_argument('--fix',
                        help=('id,value pairs to fix the flux of reaction '
                              'at certain value, multiple pairs can be '
                              'separated by ";". e.g. '
                              'id1,value1;id2,value2'))
    parser.add_argument('-o', '--output', required=True,
                        help='output file name')
    parser.add_argument('--sampler', default='optgp',
                        help='sampler, choices: [optGp (default), ACHR]')
    parser.add_argument('-n', '--nproc', type=int, default=1,
                        help=('number of processes, only works '
                              'for optGp method (default: 1)'))
    args = parser.parse_args()

    model = ModelReader.reader_from_path(args.model)
    nm = model.create_model()
    mm = nm.create_metabolic_model()

    if args.objective is not None and args.objective.lower() == 'biomass':
        args.objective = nm.biomass_reaction

    if args.sampler.lower() == 'optgp':
        args.sampler = 'optGp'
        s = optGp(mm, Solver(), fix=args.fix, objective=args.objective,
                  objective_scale=args.threshold, nproc=args.nproc)
    elif args.sampler.lower() == 'achr':
        args.sampler = 'ACHR'
        s = ACHR(mm, Solver(), fix=args.fix, objective=args.objective,
                 objective_scale=args.threshold)
    else:
        raise RuntimeError('Bad choice of sampler!')

    s.set_warmup()
    result = s.sample(args.samples)
    result.to_csv('_'.join([args.output, args.sampler, 'sampling.csv']))

    plt.figure()
    test_r = PCA(2).fit(s.warmup_flux).transform(s.warmup_flux)
    plt.scatter(test_r[:, 0], test_r[:, 1], alpha=0.1, label='warmup')
    test_r = PCA(2).fit(s.warmup_flux).transform(result)
    plt.scatter(test_r[:, 0], test_r[:, 1], alpha=0.1,
                color='r', label='sample')
    plt.legend(loc='best')
    plt.savefig('_'.join([args.output, args.sampler, 'sampling_PCA.pdf']))
