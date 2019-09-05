import argparse
import sys
from multiprocessing import Pool
from collections import OrderedDict
import pandas as pd
from scipy.linalg import null_space
import numpy as np
from psamm.lpsolver.cplex import Solver
from psamm.fluxanalysis import FluxBalanceProblem, FluxBalanceError
from psamm.datasource.native import ModelReader


class StuckError(Exception):
    def __init__(self, *args, **kwargs):
        super(StuckError, self).__init__(*args, **kwargs)


class OutOfBoundaryError(Exception):
    def __init__(self, *args, **kwargs):
        super(OutOfBoundaryError, self).__init__(*args, **kwargs)


class sampler(object):
    def __init__(self, metabolic_model, solver, fix=None,
                 objective=None, objective_scale=1, nproc=1, epsilon=1e-6):
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

    def _random_optimize(self):
        flexible = list()
        for i in self._reactions:
            limit_range = self._upper[i] - self._lower[i]
            # skip fixed reaction
            if limit_range < self._epsilon:
                continue
            flexible.append(i)
        # randomly set optimize weight
        vec = np.random.randn(len(flexible))
        # normalize to unit sphere
        vec /= np.linalg.norm(vec)
        # set up optimize object
        optimize = dict(zip(flexible, vec))
        self._p.maximize(optimize)
        xm_flux = self._get_fluxes()
        return xm_flux

    def _get_fluxes(self):
        """obtain the FBA result"""
        x = list()
        for rxn in self._reactions:
            x.append(self._p.get_flux(rxn))
        x = np.array(x)
        return x

    def set_warmup(self, more=False):
        """Set up warmup points for Monte Carlo sampling."""
        self._warmup_flux = list()
        print('Setting up warmup points...')
        # setup warmup points by maximizing and minimizing
        # the flux of each reaction
        for i in self._reactions:
            if self._upper[i] - self._lower[i] < self._epsilon:
                # skip fixed reaction
                continue
            # maximize positive flux
            try:  # maximize the flux of reaction i
                self._p.maximize(i)
                # self._p.max_min_l1(i)
                # store warmup points based on non-zero reacions only
                fluxes = self._get_fluxes()
                self._warmup_flux.append(fluxes)
            except FluxBalanceError:
                RuntimeWarning('Failed to maximize the flux of %s' % i)
                pass
            # maximize flux of reaction i in reverse direction
            try:
                self._p.maximize({i: -1})
                # self._p.max_min_l1({i: -1})
                # store warmup points based on effective reacions only
                fluxes = self._get_fluxes()
                self._warmup_flux.append(fluxes)
            except FluxBalanceError:
                RuntimeWarning('Failed to minimize the flux of %s' % i)
                pass
        if more:
            # add more warmup points by rondom optimizing
            for i in range(len(self._reactions)):
                try:
                    fluxes = self._random_optimize()
                    self._warmup_flux.append(fluxes)
                except FluxBalanceError:
                    pass
        self._warmup_flux = np.array(self._warmup_flux)
        if len(self._warmup_flux) <= 1:
            raise RuntimeError("Can't get solutions based on current "
                               "flux limitations! Please check your model.")
        print('Warmup points: %i' % len(self._warmup_flux))
        # maintain unrelated warmup points only
        u, indices = np.unique(
            np.round(self._warmup_flux, int(-np.log10(self._epsilon))),
            axis=0, return_index=True)
        # maintain un-rounded fluxes to keep accuracy
        self._warmup_flux = self._warmup_flux[sorted(indices), :]
        print(('Total reactions: %i\n'
               'Non-redundant warmup points: %i')
              % (len(self._reactions), len(self._warmup_flux)))

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
        self._name = 'optGp'

    def sample(self, nsample, k=100, maxtry=1, seed=None):
        """Artificial Centering Hit-and-Run function"""
        print('Doing artificial centering hit-and-run')
        # number of warmup points
        nwarm = len(self._warmup_flux)
        maxiter = nsample // self._nproc
        residual = nsample % self._nproc
        upper = np.array([i for i in self._upper.values()])
        lower = np.array([i for i in self._lower.values()])
        # set up random seed
        np.random.seed(seed)
        if self._nproc == 1:
            m = np.random.choice(nwarm)
            x = one_chain(m, self._warmup_flux, self._sm, self._ns, maxiter,
                          upper, lower, self._epsilon,
                          k, maxtry, seed)
            return pd.DataFrame(x, columns=self._reactions)

        elif self._nproc > 1:
            pool = Pool(self._nproc)
            start_points = np.random.choice(nwarm, self._nproc)
            tasks = task_generator(
                one_chain, start_points, self.warmup_flux, self._sm, self._ns,
                maxiter, upper, lower, self._epsilon, k, maxtry, seed,
                self._nproc, residual
            )
            x = pool.map(parallel_worker, tasks)
            pool.close()
            pool.join()
            return pd.DataFrame(np.vstack(x),
                                columns=self._reactions)


def parallel_worker(tasks):
    """Send eclapsed parameters to function"""
    func, params = tasks
    return func(*params)


def one_chain(m, warmup_flux, sm, ns, maxiter,
              upper, lower, epsilon, k, maxtry, seed=None):
    """Run one ACHR chain"""
    print('Start on point %i' % m)
    stuckcount = 0
    reprojectcount = 0
    x = list()
    nwarm = len(warmup_flux)
    # prevent altering the original warmup points
    warmup_flux = np.copy(warmup_flux)
    # get the center point
    s = warmup_flux.mean(axis=0)
    # set up starting point
    # pull back a bit to avoid stuck
    xm_flux = one_step(s, warmup_flux[m] - s, upper, lower, epsilon, 0.95)
    for niter in range(maxiter * k):
        success = False
        # wait until a successful move
        while True:
            for ntry in range(maxtry):
                n = np.random.choice(nwarm)
                xn_flux = warmup_flux[n]
                direction_flux = xn_flux - s
                try:
                    new_flux = one_step(xm_flux, direction_flux,
                                        upper, lower, epsilon)
                    success = True  # successfully got the next point
                    break
                except StuckError:
                    stuckcount += 1
                    pass
                except OutOfBoundaryError as e:
                    raise(e)
            if success:
                break
            # too many failures, move xm to new position
            # set up starting point
            # pull back to avoid stuck
            new_flux = s
            # xm_flux = one_step(s, xm_flux - s, upper, lower, epsilon, 0.9)
            # xm_flux = (xm_flux - s) * 0.9 + s
        # output only one point every k iter
        if (niter + 1) % k == 0:
            # if not np.allclose(sm.dot(new_flux), 0, 0, epsilon):
            #     reprojectcount += 1
            #     # re-project point into null space
            #     new = get_projection(new_flux, ns, epsilon, True)
            #     new_flux = ns.dot(new)
            #     # pull back a bit in case xm_flux is out of boundary
            #     while (np.sum(new_flux - upper > epsilon) > 0
            #             or np.sum(new_flux - lower < -epsilon) > 0):
            #         new_flux = one_step(xm_flux, new_flux - xm_flux,
            #                             upper, lower, epsilon, 0.95)
            # add new point
            x.append(new_flux)
            if (niter + 1) // k % 500 == 0:
                print((niter + 1) // k)
                sys.stdout.flush()
        xm_flux = new_flux
        # recalculate the center
        s = (s * (nwarm + niter) + xm_flux) / (nwarm + niter + 1)
    print('Point %i is done... %i stucks, %i reprojected into null space'
          % (m, stuckcount, reprojectcount))
    return np.array(x)


def get_projection(x, ns, epsilon, check=False):
    """obtain the FBA result in null space"""
    projection = ns.T.dot(x)
    if check:
        if np.linalg.norm(ns.dot(projection) - x, ord=2) ** 2 > epsilon:
            raise RuntimeError('Failed to project point into null space!')
        return projection
    else:
        return projection


def one_step(xm_flux, direction_flux,
             upper, lower, epsilon, step=None):
    """Move one step further on one chain"""
    if (np.sum(xm_flux - upper > epsilon) > 0
            or np.sum(xm_flux - lower < -epsilon) > 0):
        raise OutOfBoundaryError('Point is out of boundary!')
    index = np.abs(direction_flux) > epsilon
    if np.sum(index) == 0:
        raise StuckError('Direction fluxes are smaller than epsillon')
    # index = direction_flux != 0
    scale_upper = (upper
                   - xm_flux)[index] / direction_flux[index]
    scale_lower = (lower
                   - xm_flux)[index] / direction_flux[index]
    scale = np.array([scale_upper, scale_lower])
    up = scale.max(axis=0).min()
    down = scale.min(axis=0).max()
    if up - down < epsilon:
        # can't move
        raise StuckError('Cannot move!')
    # get the new point
    if step is None:
        step = np.random.uniform(down, up)
    else:
        step = down + step * (up - down)
    new_flux = xm_flux + step * direction_flux
    # got acceptable result
    return new_flux


def task_generator(one_chain, startpoints, warmup_flux, sm, ns, maxiter,
                   upper, lower, epsilon, k, maxtry, seed, limit, residual):
    for id in range(limit):
        m = startpoints[id]
        if id < residual:
            yield (one_chain,
                   (m, warmup_flux, sm, ns, maxiter + 1, upper, lower,
                    epsilon, k, maxtry, seed))
        else:
            yield (one_chain,
                   (m, warmup_flux, sm, ns, maxiter, upper, lower,
                    epsilon, k, maxtry, seed))


class allWarm(sampler):
    def __init__(self, *args, **kwargs):
        super(allWarm, self).__init__(*args, **kwargs)
        self._name = 'allWarm'

    def sample(self, nsample, k=100, maxtry=1, seed=None):
        """Artificial Centering Hit-and-Run function"""
        print('Doing artificial centering hit-and-run')
        # number of warmup points
        nwarm = len(self._warmup_flux)
        maxiter = nsample // nwarm
        residual = nsample % nwarm
        upper = np.array([i for i in self._upper.values()])
        lower = np.array([i for i in self._lower.values()])
        tasks = task_generator(
            one_chain, range(nwarm), self.warmup_flux, self._sm, self._ns,
            maxiter, upper, lower, self._epsilon, k, maxtry, seed, nwarm,
            residual)
        if self._nproc == 1:
            x = map(parallel_worker, tasks)
        elif self._nproc > 1:
            pool = Pool(self._nproc)
            x = pool.map(parallel_worker, tasks)
            pool.close()
            pool.join()
        return pd.DataFrame(np.vstack(x),
                            columns=self._reactions)


class ACHR(sampler):
    def __init__(self, *args, **kwargs):
        super(ACHR, self).__init__(*args, **kwargs)
        self._name = 'ACHR'

    def sample(self, nsample, seed=None):
        """Artificial Centering Hit-and-Run functioin"""
        print('Doing artificial centering hit-and-run')
        points_flux = [row for row in self._warmup_flux]
        npoint = len(points_flux)
        nwarm = len(self._warmup_flux)
        # center of warmup points
        s = self._warmup_flux.mean(axis=0)
        # starting point
        xm_flux = s
        upper = np.array([i for i in self._upper.values()])
        lower = np.array([i for i in self._lower.values()])
        stuckcount = 0
        reprojectcount = 0
        # set random seed
        np.random.seed(seed)
        while npoint - nwarm < nsample:
            n = np.random.choice(npoint)
            xn_flux = points_flux[n]
            direction_flux = xn_flux - s
            if np.sum(np.abs(direction_flux) > self._epsilon) == 0:
                continue
            try:
                new_flux = one_step(xm_flux, direction_flux,
                                    upper, lower, self._epsilon)
                # if not np.allclose(self._sm.dot(new_flux),
                #                    0, 0, self._epsilon):
                #     reprojectcount += 1
                #     # re-project point into null space
                #     new = get_projection(
                #         new_flux, self._ns, self._epsilon, True)
                #     new_flux = self._ns.dot(new)
                #     # pull back a bit in case xm_flux is out of boundary
                #     while (np.sum(new_flux - upper > self._epsilon) > 0
                #             or np.sum(new_flux - lower < -self._epsilon) > 0):
                #         new_flux = one_step(xm_flux, new_flux - xm_flux,
                #                             upper, lower, self._epsilon, 0.95)
                xm_flux = new_flux
                points_flux.append(xm_flux)
                # recalculate center
                s = (s * npoint + xm_flux) / (npoint + 1)
                npoint += 1
                if (npoint - nwarm) % 10000 == 0:
                    print('%i/%i' % (npoint - nwarm, nsample))
                    sys.stdout.flush()
            except StuckError:
                stuckcount += 1
                # pull back a bit
                # xm_flux = one_step(s, xm_flux - s, upper, lower,
                #                    self._epsilon, 0.9)
                # xm_flux = (xm_flux - s) * 0.9 + s
                xm_flux = s
                pass
        print('Done, %i stuck, %i projected into null space'
              % (stuckcount, reprojectcount))
        sys.stdout.flush()
        return pd.DataFrame(points_flux[nwarm:],
                            columns=self._reactions)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            'Get the distribution of reaction fluxes through '
            'Artificial Centering Hit-and-Run algorithm.'
        )
    )
    parser.add_argument('-i', '--model', help='path to metabolic model',
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
                        help=('sampler, choices: '
                              '[optGp (default), allWarm, ACHR]'))
    parser.add_argument('-n', '--nproc', type=int, default=1,
                        help=('number of processes, only works '
                              'for optGp method (default: 1)'))
    parser.add_argument('-e', '--epsilon', type=float, default=1e-6,
                        help=('precision of sampling, (default: 1e-6)'))
    parser.add_argument('--seed', type=np.int32,
                        help='random seed for sampling')
    args = parser.parse_args()

    model = ModelReader.reader_from_path(args.model)
    nm = model.create_model()
    mm = nm.create_metabolic_model()

    if args.objective is not None and args.objective.lower() == 'biomass':
        args.objective = nm.biomass_reaction

    if args.sampler.lower() == 'optgp':
        args.sampler = 'optGp'
        s = optGp(mm, Solver(), fix=args.fix, objective=args.objective,
                  objective_scale=args.threshold, nproc=args.nproc,
                  epsilon=args.epsilon)
    elif args.sampler.lower() == 'allwarm':
        args.sampler = 'allWarm'
        s = allWarm(mm, Solver(), fix=args.fix, objective=args.objective,
                    objective_scale=args.threshold, nproc=args.nproc,
                    epsilon=args.epsilon)
    elif args.sampler.lower() == 'achr':
        args.sampler = 'ACHR'
        s = ACHR(mm, Solver(), fix=args.fix, objective=args.objective,
                 objective_scale=args.threshold, epsilon=args.epsilon)
    else:
        raise RuntimeError('Bad choice of sampler!')

    s.set_warmup(more=True)
    result = s.sample(args.samples, seed=args.seed, maxtry=10)
    b = np.zeros(s._sm.shape[0])
    for index, row in result.iterrows():
        dot = s._sm.dot(row)
        if not np.allclose(dot, b, 0, 1e-6):
            print('Point %i violates equality, maximum value: %f'
                  % (index, np.max(np.abs(dot))))
    result.to_csv('_'.join([args.output, args.sampler, 'sampling.csv']))
