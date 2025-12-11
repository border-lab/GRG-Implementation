# File contains implementations for dot product on GRGs with duplicate mutations

from collections import defaultdict, deque
import numpy as np

def topo_order(grg):
    """
    Return a topological ordering of Node objects using Kahn’s algorithm.
    """
    indeg = defaultdict(int)
    children = defaultdict(list)

    for node in grg.nodes:
        for child in node.children.keys():
            indeg[child] += 1
            children[node].append(child)

    q = deque([n for n in grg.nodes if indeg[n] == 0])
    order = []

    while q:
        u = q.popleft()
        order.append(u)
        for v in children[u]:
            indeg[v] -= 1
            if indeg[v] == 0:
                q.append(v)

    if len(order) != len(grg.nodes):
        raise RuntimeError("Graph contains a cycle or invalid structure.")

    return order

def compute_descendant_leaves(grg):
    """
    Compute descendant leaf indices for every node in `grg`.

    Returns:
        leaves_by_node: dict { Node -> set of sample indices }
        hap_index: mapping { haplotype_id -> index }
        M: number of haplotypes (samples)
    """
    # Assign leaf sample indices
    hap_ids = sorted(grg.haplotype_endpoints.keys())
    hap_index = {hap_id: i for i, hap_id in enumerate(hap_ids)}
    M = len(hap_ids)

    # Topological order (children first)
    order = topo_order(grg)

    leaves_by_node = {}

    # Process in reverse topological order (roots to leaves)
    for node in reversed(order):
        # If node is a haplotype endpoint, it is a leaf sample
        if node in grg.haplotype_endpoints.values():
            # Find which haplotype ID maps to this node
            hap_id = next(h for h, n in grg.haplotype_endpoints.items() if n is node)
            leaves_by_node[node] = {hap_index[hap_id]}
        else:
            # Union of children leaf sets
            s = set()
            for child in node.children.keys():
                s |= leaves_by_node[child]
            leaves_by_node[node] = s

    return leaves_by_node, hap_index, M

def mutation_to_nodes(grg):
    """
    Build mapping mut_id -> list of Node objects where it appears.
    """
    m2nodes = defaultdict(list)
    for node in grg.nodes:
        for m in node.mut:
            if m is not None:
                m2nodes[m].append(node)
    return m2nodes

def xu_with_duplicates(grg, u):
    """
    Compute Xu where:
      - u has length N_mutations (N == number of unique mutation IDs in GRG)
      - each mutation may appear multiple times in the GRG

    Returns:
        x_u: numpy array of length M (haplotype/sample vector)
    """
    leaves_by_node, hap_index, M = compute_descendant_leaves(grg)
    m2nodes = mutation_to_nodes(grg)

    # Build stable ordering of mutation IDs and mapping to indices of u
    mut_ids = sorted(m2nodes.keys())
    if len(mut_ids) != len(u):
        raise ValueError(
            f"Length of u ({len(u)}) does not match the number of mutation IDs ({len(mut_ids)})."
        )
    mut_id_to_u_index = {m: i for i, m in enumerate(mut_ids)}

    x_u = np.zeros(M, dtype=float)

    for mut_id, nodes in m2nodes.items():
        idx = mut_id_to_u_index[mut_id]
        uj = u[idx]
        if uj == 0:
            continue
        union_leaves = set()
        for node in nodes:
            union_leaves |= leaves_by_node[node]
        for leaf in union_leaves:
            x_u[leaf] += uj

    return x_u


def xtv_with_duplicates(grg, v):
    """
    Compute X^T v where:
      - v has length M haplotypes/samples
      - mutations may appear multiple times

    Returns:
        (xtv, mut_ids): xtv is numpy array of length N (same order as mut_ids)
                         mut_ids is the list of mutation IDs corresponding to xtv entries
    """
    leaves_by_node, hap_index, M = compute_descendant_leaves(grg)
    m2nodes = mutation_to_nodes(grg)

    if len(v) != M:
        raise ValueError(
            f"Length of v ({len(v)}) does not match number of haplotypes M ({M})."
        )

    mut_ids = sorted(m2nodes.keys())
    xtv = np.zeros(len(mut_ids), dtype=float)

    for out_idx, mut_id in enumerate(mut_ids):
        nodes = m2nodes[mut_id]
        union_leaves = set()
        for node in nodes:
            union_leaves |= leaves_by_node[node]
        total = sum(v[leaf] for leaf in union_leaves)
        xtv[out_idx] = total

    return xtv