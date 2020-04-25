#!/usr/bin/env python
# coding: utf-8
"""Parse desired reaction and metabolite information from BiGG Models.

Using API provided by BiGG Models (http://bigg.ucsd.edu/data_access).
get_metabolite(id, model, delay) for parsing metabolite
get_reaction(id, model, delay) for parsing reaction
"""

import json
import urllib
from time import sleep


def parse_url(id, model, type):
    """Parse url for query id, return dictionary.

    The value for type is either 'metabolites' or 'reactions'.
    """
    apiurl = 'http://bigg.ucsd.edu/api/v2'
    if model == 'universal':
        url = '/'.join([apiurl, model, type, id])
    else:
        url = '/'.join([apiurl, 'models', model, type, id])
    req = urllib.request.urlopen(url)
    return json.load(req)


def get_metabolite(id, model='universal', delay=0.1):
    """Get query metabolite."""
    sleep(delay)  # set a delay to not exceed 10 requres per second
    result = parse_url(id, model, 'metabolites')
    if model == 'universal':
        formula = result['formulae'][0]
        charge = result['charges'][0]
    else:
        formula = result['formula']
        charge = result['charge']
    return {'id': result['bigg_id'],
            'name': result['name'],
            'formula': formula,
            'charge': charge,
            'kegg': [i['id']
                     for i in result['database_links']['KEGG Compound']],
            'biocyc': [i['id'] for i in result['database_links']['BioCyc']]}
