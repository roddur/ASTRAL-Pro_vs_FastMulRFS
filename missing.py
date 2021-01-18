'''
Written by Mrinmoy S. Roddur on December, 2020

Deletes clades from gene trees based on M_clade model described in 
Nute, Michael, et al. "The performance of coalescent-based species tree estimation methods under models of missing data."

Description: 
First, it lists all the clades from the species tree having the number of species above and below 
two thresholds. Next, the gene trees from which a clade cannot be deleted (have enough leaves and a clade 
has been found contained entirely in the gene tree after a specified number of tries) are filtered out.
From the remaining gene trees, a certain percentage of all trees have been chosen and from each tree, a random clade listed from the 
species tree has been deleted. The process is the same as described as Nute et. al. when all the 
gene trees have the same number of species.

'''
import treeswift
import random
import argparse
import numpy as np

def clade(node):
    leaves = set()
    for i in node.traverse_leaves():
        leaves.add(i.get_label().split('_')[0])
    return frozenset(leaves)

def main(args):

    s_tree = treeswift.read_tree_newick(args.speciestree)
    total_taxa = len(clade(s_tree.root))
    all_clades = set()

    for i in s_tree.traverse_preorder(leaves=False):
        for j in i.child_nodes():
            newclade = clade(j)
            if total_taxa*args.lower < len(newclade) < total_taxa*args.upper:
                all_clades.add(newclade)
    
    g_trees = treeswift.read_tree_newick(args.genetree)
    for g_tree in g_trees:
        if len(clade(g_tree.root)) > total_taxa*args.lower and random.random() < args.missing:
            for try_ in range(args.tries):
                gtreebak = treeswift.read_tree_newick(g_tree.newick())
                selected_clade = random.sample(all_clades, 1)[0]
                for i in g_tree.traverse_preorder(leaves=False):
                    for j in i.child_nodes():
                        if len(selected_clade.intersection(clade(j))) == 0:
                            i.remove_child(j)
                    if i.num_children() == 1:
                        i.contract()
                if g_tree.root.num_children() == 0 or len(clade(g_tree.root)) < total_taxa*args.lower:
                    g_tree = gtreebak
                else:
                    break
        print(g_tree.newick())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-g", "--genetree", type=str, help="Input gene trees", required=True)
    parser.add_argument("-s", "--speciestree", type=str, help="Input species tree", required=True)
    parser.add_argument("-l", "--lower", type=float, help="Lower limit of the percentege of taxon (excluding)", default=0.2)
    parser.add_argument("-u", "--upper", type=float, help="Upper limit of the percentege of taxon (excluding)", default=1)
    parser.add_argument("-m", "--missing", type=float, help="Percentage of trees with missing taxa", default=0.95)
    parser.add_argument("-t", "--tries", type=float, help="Number of tries to get a gene tree with missing data", default=5)
    main(parser.parse_args())
