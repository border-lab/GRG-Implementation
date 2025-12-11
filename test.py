import json

import dot_product
import grg
import numpy as np
from recombination import *
from plot_grg import *
from dot_product_dup import xu_with_duplicates, xtv_with_duplicates

n0 = grg.Node(0, (2, 3))
n1 = grg.Node(1, (1))
n2 = grg.Node(2, (4))
n3 = grg.Node(3, (5))
n4 = grg.Node(4, (6, 7))
n5 = grg.Node(5, (8))
n6 = grg.Node(6, (9))
n7 = grg.Node(7, (10, 11))

full = [0, 12]

n6.connect(n0, full)
n6.connect(n5, full)
n7.connect(n0, full)
n7.connect(n4, full)
n5.connect(n1, full)
n5.connect(n4, full)
n4.connect(n2, full)
n4.connect(n3, full)

all_nodes = {n0, n1, n2, n3, n4, n5, n6, n7}

endpoints_map = {
    "hapA": n0,  # Path: 6 -> 0 or 7 -> 0
    "hapB": n1,  # Path: 6 -> 5 -> 1
    "hapC": n2,  # Path: 6 -> 5 -> 4 -> 2
    "hapD": n3,  # Path: 7 -> 4 -> 3
}
individual_map = {
    "hapA": "Individual_1",
    "hapB": "Individual_1",
    "hapC": "Individual_2",
    "hapD": "Individual_2",
}

grg_ = grg.GRG(all_nodes, endpoints_map, individual_map)
# plot_grg(grg_)
mutation_map = grg_.get_mutation_to_individuals()

individual_to_mutation_map = grg_.get_individuals_to_mutations()
print(json.dumps(individual_to_mutation_map, indent=2))

# plot_grg(grg_)
print("=== BEFORE RECOMBINATION ===")
print(json.dumps(mutation_map, indent=2))
print(json.dumps(individual_to_mutation_map, indent=2))

# Recombine hapA and hapC at breakpoint 6 into new hapE
offspring_node = recombine(
    grg_,
    hapA_id="hapB",
    hapB_id="hapA",
    breakpoint=3,
    new_hap_id="hapE",
    new_ind_id="Individual_3",
    ADDING_NODES=True,
)

print("\nNew offspring haplotype 'hapE' is node:", offspring_node)

# Recompute maps after recombination
mutation_map_after = grg_.get_mutation_to_individuals()
individual_map_after = grg_.get_individuals_to_mutations()

# Dot Product Tests
X = grg_.to_matrix()
u = np.array([1] * X.shape[1], dtype=float)  # Vector of ones
v = np.array([1] * X.shape[0], dtype=float)  # Vector of ones

print("\n=== DOT PRODUCT TESTS ===")
print(np.dot(X, u))  
print(np.dot(X.T, v))

# Results
print("\n=== DOT PRODUCT WITH DUPLICATES TESTS ===")
print("X u with duplicates:", xu_with_duplicates(grg_, u))
print("X^T v with duplicates:", xtv_with_duplicates(grg_, v))


print("\n=== AFTER RECOMBINATION ===")
print(json.dumps(mutation_map_after, indent=2))
print(json.dumps(individual_map_after, indent=2))
plot_grg(grg_)
