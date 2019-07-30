# coding: utf-8
from psamm.datasource.native import ModelReader
from psamm.fluxanalysis import FluxBalanceProblem
from psamm.lpsolver.cplex import Solver
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--model', help='Model file to be loaded.', required=True)
parser.add_argument('-f', '--fixed-reaction', action='append',
                    help=('"reaction,value" pair to fix the flux of reaction. '
                          'Can be applied multiple times to fix multiple '
                          'reactions.'))
parser.add_argument('-r', '--varying-reaction',
                    help='Reaction whose flux should be varied.',
                    required=True)
parser.add_argument('--min', type=float,
                    help='Min flux for the varying reaction.')
parser.add_argument('--max', type=float,
                    help='Max flux for the varying reaction.')
parser.add_argument('-s', '--steps', type=int,
                    help=('How many steps are between min '
                          'and max flux of the varying'
                          'reaction (default = 10).'),
                    default=10)
parser.add_argument('--objective', help='Reaction to use as objective.')

args = parser.parse_args()

model = ModelReader.reader_from_path(args.model)
nm = model.create_model()
mm = nm.create_metabolic_model()

p = FluxBalanceProblem(mm, Solver())

if not mm.has_reaction(args.varying_reaction):
    sys.exit('Error: %s is not in model' % (args.varying_reaction))
if not mm.has_reaction(args.objective):
    sys.exit('Error: %s is not in model' % (args.objective))

fvs = list()
for i in args.fixed_reaction:
    r, v = i.split(',')
    if not mm.has_reaction(r):
        sys.exit('Error: %s is not in model' % (r))
    fvs.append(p.get_flux_var(r) == float(v))

# set fixed reactions as linear constraints
constraints = p.prob.add_linear_constraints(*fvs)

if args.max is None:
    p.maximize(args.varying_reaction)
    flux_max = p.get_flux(args.varying_reaction)
else:
    flux_max = args.max

if args.min is None:
    p.maximize({args.varying_reaction: -1})
    flux_min = p.get_flux(args.varying_reaction)
else:
    flux_min = args.min

# varying reaction to get the result
for i in range(args.steps):
    fixed_flux = flux_min + i * (flux_max - flux_min) / float(args.steps - 1)
    fvs.append(p.get_flux_var(args.varying_reaction) == fixed_flux)
    constraints = p.prob.add_linear_constraints(*fvs)
    p.maximize(args.objective)
    # print results for one iter
    for rxn in mm.reactions:
        value = p.get_flux(rxn)
        print('%s\t%f\t%.5f' % (rxn, fixed_flux, value))
    # remove constraint
    fvs.pop()
    c = constraints.pop()
    c.delete()
