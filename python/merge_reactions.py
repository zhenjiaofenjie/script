from psamm.datasource.native import (parse_reaction_file as pr,
                                     ModelWriter)
import argparse


def merge(template, addition):
    template_rxn = pr(template)
    addition_rxn = pr(addition)
    union = list(template_rxn)
    ids = [r.id for r in union]
    for r in addition_rxn:
        if r.id not in ids:
            union.append(r)
    return union


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=('Add novel reactions from "additional" file to '
                     '"template" file, ignore reactions based on reaction id.')
    )
    parser.add_argument('--template', help='template reaction YAML file')
    parser.add_argument('--addition', help='additional reaction YAML file')
    parser.add_argument('--out', help='output YAML file')
    args = parser.parse_args()

    union = merge(args.template, args.addition)
    with open(args.out, 'w') as o:
        ModelWriter().write_reactions(o, union)
