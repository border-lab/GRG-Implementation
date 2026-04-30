"""
grg_numpy_baseline.py

Contains the baseline approach for extracting the GRG into a dense NumPy matrix.
Uses a robust bottom-up traversal to bypass PyGRGL's internal root-caching bugs.
"""

import numpy as np
from collections import deque

def estimate_numpy_memory(num_mutations, num_samples, dtype_size=1):
    """Estimates the memory required to instantiate the dense NumPy matrix in MB."""
    bytes_required = num_mutations * num_samples * dtype_size
    return bytes_required / (1024 * 1024)

def grg_to_numpy(g):
    """
    Converts a GRG to a dense NumPy genotype matrix using a bottom-up traversal.
    """
    # 1. Map samples to columns
    sample_nodes = [i for i in range(g.num_nodes) if g.is_sample(i)]
    sample_nodes.sort()
    node_to_col = {node: idx for idx, node in enumerate(sample_nodes)}
    
    # 2. Bottom-up traversal: For every node, which samples descend from it?
    node_samples = {i: set() for i in range(g.num_nodes)}
    
    for sample in sample_nodes:
        queue = deque([sample])      # use deque
        visited = {sample}
        
        while queue:
            curr = queue.popleft()   # O(1) instead of pop(0)
            node_samples[curr].add(sample)
            
            # Walk UP the graph to all ancestors
            for parent in g.get_up_edges(curr):
                if parent not in visited:
                    visited.add(parent)
                    queue.append(parent)
    
    # 3. Initialize dense matrix
    genotype_matrix = np.zeros((g.num_mutations, g.num_samples), dtype=np.int8)
    
    # 4. Populate matrix
    for node_id, mutation_id in g.get_node_mutation_pairs():
        # Only populate if the node actually has descendant samples
        if node_id in node_samples:
            for sample_node in node_samples[node_id]:
                genotype_matrix[mutation_id, node_to_col[sample_node]] = 1
            
    return genotype_matrix.T


def grg_to_numpy_parallel(g, n_jobs=4):
    """
    Parallelized version of grg_to_numpy using joblib.
    """
    try:
        from joblib import Parallel, delayed
    except ImportError:
        return grg_to_numpy(g)
    
    sample_nodes = [i for i in range(g.num_nodes) if g.is_sample(i)]
    sample_nodes.sort()
    node_to_col = {node: idx for idx, node in enumerate(sample_nodes)}
    
    # Bottom-up traversal
    node_samples = {i: set() for i in range(g.num_nodes)}
    for sample in sample_nodes:
        queue = [sample]
        visited = set([sample])
        while queue:
            curr = queue.pop(0)
            node_samples[curr].add(sample)
            for parent in g.get_up_edges(curr):
                if parent not in visited:
                    visited.add(parent)
                    queue.append(parent)
    
    mutation_pairs = list(g.get_node_mutation_pairs())
    
    def process_chunk(chunk):
        rows, cols = [], []
        for node_id, mutation_id in chunk:
            if node_id in node_samples:
                for sample_node in node_samples[node_id]:
                    rows.append(mutation_id)
                    cols.append(node_to_col[sample_node])
        return rows, cols
    
    chunk_size = max(1, len(mutation_pairs) // n_jobs)
    chunks = [mutation_pairs[i:i+chunk_size] 
              for i in range(0, len(mutation_pairs), chunk_size)]
    
    results = Parallel(n_jobs=n_jobs)(delayed(process_chunk)(chunk) for chunk in chunks)
    
    rows = np.concatenate([r[0] for r in results]) if results else np.array([], dtype=int)
    cols = np.concatenate([r[1] for r in results]) if results else np.array([], dtype=int)
    
    genotype_matrix = np.zeros((g.num_mutations, g.num_samples), dtype=np.int8)
    if len(rows) > 0:
        genotype_matrix[rows, cols] = 1
    
    return genotype_matrix.T

def simulate_numpy_recombination(benchmark, genotype_matrix, bp_range, expected_crossovers=1.5):
    """
    Performs true recombination using standard NumPy array slicing.
    Each parent mates exactly once, and each pair produces 2 children.
    """
    num_samples, genome_length = genotype_matrix.shape
    total_bp = 0
    
    # Calculate exact offspring count (handles odd numbers safely by flooring pairs)
    num_pairs = num_samples // 2
    num_offspring = num_pairs * 2
    # print(f"Simulating recombination for {num_samples} samples ({num_pairs} pairs) to produce {num_offspring} offspring.")
    
    offspring_matrix = np.zeros((num_offspring, genome_length), dtype=genotype_matrix.dtype)
    
    # Shuffle parent indices for random pairing
    parent_indices = np.arange(num_samples)
    np.random.shuffle(parent_indices)
    
    offspring_idx = 0
    
    # Iterate through pairs
    for i in range(0, num_samples - 1, 2):
        p1 = parent_indices[i]
        p2 = parent_indices[i+1]
        
        # Each pair produces 2 children
        for j in range(2):
            if offspring_idx >= num_offspring:
                return offspring_matrix, total_bp
            bp, num_bp = benchmark.get_breakpoints(bp_range, expected_crossovers)
            total_bp += num_bp
            # print("testcase breakpoints:", offspring_idx, j)  # Debug print for breakpoints
            # bp = testcases[offspring_idx][j] # Use pre-generated breakpoints for consistency
            current_parent = p1 if np.random.random() < 0.5 else p2
            other_parent = p2 if current_parent == p1 else p1
            
            start_idx = 0
            for b in bp:
                offspring_matrix[offspring_idx, start_idx:b] = genotype_matrix[current_parent, start_idx:b]
                start_idx = b
                current_parent, other_parent = other_parent, current_parent
                
            offspring_matrix[offspring_idx, start_idx:] = genotype_matrix[current_parent, start_idx:]
            offspring_idx += 1
        
        if offspring_idx >= num_offspring:
            return offspring_matrix, total_bp
            
    return offspring_matrix, total_bp