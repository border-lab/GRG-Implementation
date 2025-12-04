import math
from grg import *

def get_overlap(a, b):
    l = max(a[0], b[0])
    r = min(a[1], b[1])
    return [l, r] if l < r else None

def get_relevant_mutations(muts, L, R):
    return sorted([m for m in muts if m is not None and L <= m < R])

def get_mutation_bounds_above(node):
    """
    Get (min_mut, max_mut) for all mutations reachable above this node.
    Traverses all ancestors and returns the min and max mutation numbers.
    """
    visited = set()
    all_muts = []
    stack = [node]
    
    while stack:
        current = stack.pop()
        if current in visited:
            continue
        visited.add(current)
        
        # Collect mutations from current node
        all_muts.extend([m for m in current.mut if m is not None])
        
        # Add all parents to stack
        for parent in current.parents.keys():
            if parent not in visited:
                stack.append(parent)
    
    if not all_muts:
        return None
    return (min(all_muts), max(all_muts))


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
    
    print("parent edges:", parent_edges)
    
    if not parent_edges:
        current_node.connect(offspring, [L, R])
        return

    for parent, overlap in parent_edges:
        recurse_attach(grg, parent, offspring, overlap, adding_nodes)


def recurse_duplication(grg, current_node, offspring, interval):
    L, R = interval
    if L >= R:
        return
    
    # Get mutations from current node that fall within [L, R]
    relevant_muts = get_relevant_mutations(current_node.mut, L, R)
    
    # Get mutation bounds reachable above this node
    bounds_above = get_mutation_bounds_above(current_node)
    
    # Condition 1: Both current node mutations and ancestors' mutations are within [L, R]
    if relevant_muts:
        current_mut_set = set(current_node.mut)
        relevant_mut_set = set(relevant_muts)
        
        if current_mut_set == relevant_mut_set:  # All current mutations are needed
            if bounds_above is not None:
                min_above, max_above = bounds_above
                if min_above >= L and max_above < R:  # All ancestors' mutations are within [L, R]
                    # Directly connect current node to offspring
                    current_node.connect(offspring, [L, R])
                    return
    
    # Condition 2: Mutations above contain mutations we don't need
    if bounds_above is not None:
        min_above, max_above = bounds_above
        # Check if bounds_above are NOT within [L, R]
        if not (min_above >= L and max_above < R):
            # Bounds above are not a subset of [L, R], so we need to duplicate
            # Copy over mutations from current node that are within [L, R]
            if relevant_muts:
                for m in relevant_muts:
                    offspring.mut = tuple(sorted(offspring.mut + (m,)))
            
            # Recursively process parent nodes
            for parent in current_node.parents.keys():
                recurse_duplication(grg, parent, offspring, [L, R])
            return
    
    # Condition 3: Bounds above are within [L, R], but not all current mutations are needed
    if bounds_above is not None:
        min_above, max_above = bounds_above
        if min_above >= L and max_above < R:  # Bounds above are within [L, R]
            current_mut_set = set(current_node.mut)
            relevant_mut_set = set(relevant_muts) if relevant_muts else set()
            
            if current_mut_set != relevant_mut_set and relevant_mut_set:
                # Not all current mutations are needed, copy needed ones
                for m in relevant_muts:
                    offspring.mut = tuple(sorted(offspring.mut + (m,)))
                
                # Attach offspring to all parents with original interval
                for parent in current_node.parents.keys():
                    offspring.connect(parent, [L, R])
                return
    
    # Handle nodes with no parents
    if not current_node.parents:
        current_mut_set = set(current_node.mut)
        relevant_mut_set = set(relevant_muts) if relevant_muts else set()
        
        if current_mut_set == relevant_mut_set and relevant_mut_set:
            # All mutations are within [L, R], connect directly
            current_node.connect(offspring, [L, R])
        else:
            # Not all mutations are within [L, R], copy the needed ones
            if relevant_muts:
                for m in relevant_muts:
                    offspring.mut = tuple(sorted(offspring.mut + (m,)))
        return
    
    # Default: continue traversing up to parents
    for parent in current_node.parents.keys():
        recurse_duplication(grg, parent, offspring, [L, R])


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