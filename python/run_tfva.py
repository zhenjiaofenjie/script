#!/usr/bin/env python
# coding: utf-8

# In[3]:


from psamm.datasource.native import ModelReader
from psamm.fluxanalysis import FluxBalanceProblem
from psamm.lpsolver.cplex import Solver
# from psamm.fluxanalysis import flux_variability
import sys
import timeout_decorator


# In[4]:


# rxn_list = []
# with open(sys.argv[2]) as f:
#     for row in f:
#         if '#' in row:
#             continue
#         rxn_list.append(row.strip())


# In[5]:


model = ModelReader.reader_from_path(sys.argv[1])
nm = model.create_model()
mm = nm.create_metabolic_model()

p = FluxBalanceProblem(mm, Solver())
p.maximize(nm.biomass_reaction)
threshold = p.get_flux(nm.biomass_reaction)
p.add_thermodynamic()
p.prob.add_linear_constraints(
    p.get_flux_var(nm.biomass_reaction) >= threshold)

rxn_list = [rxn for rxn in mm.reactions]
# In[ ]:


@timeout_decorator.timeout(10, use_signals=False)
def get_bound(p, rxn):
    p.maximize(rxn)
    upper = p.get_flux(rxn)
    p.maximize({rxn: -1})
    lower = p.get_flux(rxn)
    return lower, upper


# In[ ]:
total = len(rxn_list)
count = 0
with open(sys.argv[2] + '.out', 'w') as o:
    for rxn in rxn_list:
        count += 1
        print('Solving %i/%i...' % (count, total))
        try:
            lower, upper = get_bound(p, rxn)
        except timeout_decorator.TimeoutError:
            print("Can't solve %s" % rxn)
            continue
        o.write('%s\t%f\t%f\n' % (rxn, lower, upper))
