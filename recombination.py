import math
from grg import *

def get_overlap(a, b):
    l = max(a[0], b[0])
    r = min(a[1], b[1])
    return [l, r] if l < r else None

def get_relevant_mutations(muts, L, R):
    return sorted([m for m in muts if m is not None and L <= m < R])

def insert_middle_node(grg, original_node, muts_to_extract):
    current_max = max((n.n for n in grg.nodes if n.n >= 0), default=0)
    new_idx = current_max + 1
    
    middle_node = Node(new_idx, tuple(muts_to_extract))
    grg.nodes.add(middle_node)

    existing_parents = list(original_node.parents.items())
    
    for parent, intervals in existing_parents:
        for iv in intervals:
            parent.connect(middle_node, iv)
            if original_node in parent.children:
                del parent.children[original_node]
        
    original_node.parents = {} 
    middle_node.connect(original_node, [0, math.inf])

    current_set = set(original_node.mut)
    remove_set = set(muts_to_extract)
    original_node.mut = tuple(sorted(list(current_set - remove_set)))

    return middle_node

def recurse_attach(grg, current_node, offspring, interval, adding_nodes):
    L, R = interval
    if L >= R: return

    relevant_muts = get_relevant_mutations(current_node.mut, L, R)

    if relevant_muts:
        target_node = current_node
        
        all_muts = set(current_node.mut)
        rel_muts_set = set(relevant_muts)
        
        if adding_nodes and (all_muts != rel_muts_set):
            target_node = insert_middle_node(grg, current_node, relevant_muts)
        
        target_node.connect(offspring, [L, R])
        return
    
    parent_edges = []
    for parent, intervals in current_node.parents.items():
        for iv in intervals:
            overlap = get_overlap(iv, [L, R])
            if overlap:
                parent_edges.append((parent, overlap))
    
    if not parent_edges:
        current_node.connect(offspring, [L, R])
        return

    for parent, overlap in parent_edges:
        recurse_attach(grg, parent, offspring, overlap, adding_nodes)

def recombine(grg_obj, hapA_id, hapB_id, breakpoint, new_hap_id, new_ind_id, ADDING_NODES=True):
    current_min = min((n.n for n in grg_obj.nodes), default=0)
    new_idx = min(current_min, 0) - 1
    
    offspring = Node(new_idx, ()) 
    grg_obj.add_node(offspring)
    
    grg_obj.haplotype_endpoints[new_hap_id] = offspring
    grg_obj.haplotype_to_individual_map[new_hap_id] = new_ind_id
    
    max_pos = 0
    for n in grg_obj.nodes:
        if n.mut and n.mut != (None,):
            max_pos = max(max_pos, max(n.mut))
    max_limit = max_pos + 100

    nodeA = grg_obj.haplotype_endpoints[hapA_id]
    nodeB = grg_obj.haplotype_endpoints[hapB_id]

    if breakpoint > 0:
        recurse_attach(grg_obj, nodeA, offspring, [0, breakpoint], ADDING_NODES)
    
    if breakpoint < max_limit:
        recurse_attach(grg_obj, nodeB, offspring, [breakpoint, max_limit], ADDING_NODES)

    return offspring