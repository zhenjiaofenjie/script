import argparse
import sys
from multiprocessing import Pool
from collections import OrderedDict
import pandas as pd
from scipy.linalg import null_space, lstsq
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
                 objective=None, objective_scale=1, nproc=1, epsilon=1e-9):
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
        if check:
            projection = lstsq(self._ns, x)
            if projection[1] > self._epsilon:
                raise RuntimeError('Point violate null space!')
        else:
            projection = self._ns.T.dot(x)
        return projection[0]

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
        return xm_flux

    def _get_fluxes(self):
        """obtain the FBA result"""
        x = list()
        for rxn in self._reactions:
            x.append(self._p.get_flux(rxn))
        x = np.array(x)
        return x

    def set_warmup(self):
        """Set up warmup points for Monte Carlo sampling."""
        self._warmup_flux = list()
        print('Setting up warmup points...')
        for i in self._reactions:
            if self._upper[i] - self._lower[i] < self._epsilon:
                # skip fixed reaction
                continue
            # maximize positive flux
            if self._upper[i] > self._epsilon:
                try:  # maximize the flux of reaction i
                    self._p.maximize(i)
                    # store warmup points based on non-zero reacions only
                    fluxes = self._get_fluxes()
                    self._warmup_flux.append(fluxes)
                except FluxBalanceError:
                    pass
            # maximize flux of reaction i in reverse direction
            if self._lower[i] < -self._epsilon:
                try:
                    self._p.maximize({i: -1})
                    # store warmup points based on effective reacions only
                    fluxes = self._get_fluxes()
                    self._warmup_flux.append(fluxes)
                except FluxBalanceError:
                    pass
        self._warmup_flux = np.array(self._warmup_flux)
        if len(self._warmup_flux) <= 1:
            raise RuntimeError("Can't get solutions based on current "
                               "flux limitations! Please check your model.")
        print('Warmup points: %i' % len(self._warmup_flux))
        # maintain unrelated warmup points only
        self._warmup_flux = np.unique(self._warmup_flux, axis=0)
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

    def sample(self, nsample, k=100, maxtry=1):
        """Artificial Centering Hit-and-Run functioin"""
        print('Doing artificial centering hit-and-run')
        # number of warmup points
        nwarm = len(self._warmup_flux)
        maxiter = nsample * k // self._nproc + 1
        upper = np.array([i for i in self._upper.values()])
        lower = np.array([i for i in self._lower.values()])
        if self._nproc == 1:
            m = np.random.choice(nwarm)
            x = one_chain(m, self._warmup_flux, self._ns, maxiter,
                          upper, lower, self._epsilon,
                          k, maxtry)
            return pd.DataFrame(x, columns=self._reactions)

        elif self._nproc > 1:
            pool = Pool(self._nproc)
            tasks = ((one_chain,
                      (m, self._warmup_flux, self._ns, maxiter,
                       upper, lower, self._epsilon, k, maxtry))
                     for m in np.random.choice(nwarm, self._nproc))
            x = pool.map(parallel_worker, tasks)
            pool.close()
            pool.join()
            return pd.DataFrame(np.vstack(x),
                                columns=self._reactions)


def parallel_worker(tasks):
    """Send eclapsed parameters to function"""
    func, params = tasks
    return func(*params)


def one_chain(m, warmup_flux, ns, maxiter, upper, lower, epsilon, k, maxtry):
    """Run one ACHR chain"""
    print('Start on point %i' % m)
    x = list()
    nwarm = len(warmup_flux)
    # prevent altering the original warmup points
    warmup_flux = np.copy(warmup_flux)
    # get the center point
    s = warmup_flux.mean(axis=0)
    # set up starting point
    # pull back a bit to avoid stuck
    xm_flux = one_step(s, warmup_flux[m] - s, upper, lower, epsilon, 0.9)
    for niter in range(maxiter):
        success = False
        # wait until a successful move
        while True:
            for ntry in range(maxtry):
                n = np.random.choice(nwarm)
                xn_flux = warmup_flux[n]
                direction_flux = xn_flux - s
                try:
                    xm_flux = one_step(xm_flux, direction_flux,
                                       upper, lower, epsilon)
                    success = True  # successfully got the next point
                    break
                except StuckError:
                    pass
                except OutOfBoundaryError:
                    pass
            if success:
                break
            # too many failures, move xm to new position
            # set up starting point
            # pull back to avoid stuck
            xm_flux = s
            # xm_flux = one_step(s, xm_flux - s, upper, lower, epsilon, 0.9)
            # xm_flux = (xm_flux - s) * 0.9 + s
        # recalculate the center
        s = (s * (nwarm + niter) + xm_flux) / (nwarm + niter + 1)
        # output only one point every k iter
        if (niter + 1) % k == 0:
            # if not np.allclose(self._sm.dot(xm_flux), 0, 0, epsilon):
            #     # re-project point into null space
            #     xm = get_projection(xm_flux, ns, epsilon)
            #     xm_flux = ns.dot(xm)
            # add new point
            x.append(xm_flux)
            if (niter + 1) // k % 500 == 0:
                print((niter + 1) // k)
                sys.stdout.flush()
    print('Point %i is done...' % m)
    return np.array(x)


def get_projection(x, ns, epsilon):
    """obtain the FBA result in null space"""
    projection = lstsq(ns, x)
    if projection[1] > epsilon:
        raise RuntimeError('Failed to project point into null space!')
    return projection[0]


def one_step(xm_flux, direction_flux,
             upper, lower, epsilon, step=None):
    """Move one step further on one chain"""
    if (np.sum(xm_flux - upper > epsilon) > 0
            or np.sum(xm_flux - lower < -epsilon) > 0):
        raise OutOfBoundaryError('Point is out of boundary!')
    index = np.abs(direction_flux) > epsilon
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


class ACHR(sampler):
    def __init__(self, *args, **kwargs):
        super(ACHR, self).__init__(*args, **kwargs)

    def sample(self, nsample):
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
        while npoint - nwarm < nsample:
            n = np.random.choice(npoint)
            xn_flux = points_flux[n]
            direction_flux = xn_flux - s
            try:
                xm_flux = one_step(xm_flux, direction_flux,
                                   upper, lower, self._epsilon)
                # # reproject
                # xm = self._ns.T.dot(xm_flux)
                # xm_flux = self._ns.dot(xm)
                points_flux.append(xm_flux)
                # recalculate center
                s = (s * npoint + xm_flux) / (npoint + 1)
                npoint += 1
                if (npoint - nwarm) % 10000 == 0:
                    print('%i/%i' % (npoint - nwarm, nsample))
                    sys.stdout.flush()
            except StuckError:
                # pull back a bit
                # xm_flux = one_step(s, xm_flux - s, upper, lower,
                #                    self._epsilon, 0.9)
                # xm_flux = (xm_flux - s) * 0.9 + s
                xm_flux = s
                pass
            except OutOfBoundaryError:
                xm_flux = s
        return pd.DataFrame(points_flux[nwarm:],
                            columns=self._reactions)


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
    parser.add_argument('-e', '--epsilon', type=float, default=1e-9,
                        help=('precision of sampling, (default: 1e-9)'))
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
    elif args.sampler.lower() == 'achr':
        args.sampler = 'ACHR'
        s = ACHR(mm, Solver(), fix=args.fix, objective=args.objective,
                 objective_scale=args.threshold, epsilon=args.epsilon)
    else:
        raise RuntimeError('Bad choice of sampler!')

    s.set_warmup()
    result = s.sample(args.samples)
    b = np.zeros(s._sm.shape[0])
    for index, row in result.iterrows():
        dot = s._sm.dot(row)
        if not np.allclose(dot, b, 0, 1e-5):
            print('Point %i violates equality, maximum value: %f'
                  % (index, np.max(np.abs(dot))))
    result.to_csv('_'.join([args.output, args.sampler, 'sampling.csv']))
