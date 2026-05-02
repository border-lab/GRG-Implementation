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


def _fill_one_offspring_alternating(matrix, p1, p2, bp, dest_row, dest, shuffle_start_parent):
    """One offspring: shared ``bp``; alternate p1/p2 along segments (random start if True)."""
    if shuffle_start_parent:
        current = p1 if np.random.random() < 0.5 else p2
    else:
        current = p1
    other = p2 if current == p1 else p1
    start_idx = 0
    for b in bp:
        dest[dest_row, start_idx:b] = matrix[current, start_idx:b]
        start_idx = b
        current, other = other, current
    dest[dest_row, start_idx:] = matrix[current, start_idx:]


def _fill_shared_breakpoints_k_offspring(matrix, p1, p2, bp, k, row_base, dest):
    """
    One shared breakpoint set per couple. For each segment, sibling rows partition the
    two parents: slot ``(j + segment_index + phase) % k`` maps to p1 if
    ``slot < k // 2`` else p2, with random ``phase`` in ``[0, k)`` per couple.

    For k==2 this is complementary offspring (no duplicate parental segment assignment);
    for larger k, each segment splits offspring between p1 and p2 as evenly as
    possible (floor(k/2) vs ceil(k/2)), rotated along the chromosome.

    k==1 uses alternating segments along the shared breakpoints only (single gamete).
    """
    genome_length = matrix.shape[1]
    bp = np.asarray(bp, dtype=int)

    if k == 1:
        _fill_one_offspring_alternating(
            matrix, p1, p2, bp, row_base, dest, shuffle_start_parent=True
        )
        return

    phase = np.random.randint(0, k)
    half = k // 2
    segment_starts = [0] + [int(x) for x in bp]
    segment_ends = list(bp) + [genome_length]
    num_segments = len(segment_starts)

    for s in range(num_segments):
        a, bnd = segment_starts[s], segment_ends[s]
        if a >= bnd:
            continue
        for j in range(k):
            slot = (j + s + phase) % k
            parent = p1 if slot < half else p2
            dest[row_base + j, a:bnd] = matrix[parent, a:bnd]


def simulate_numpy_recombination(benchmark, genotype_matrix, bp_range, expected_crossovers=1.5):
    """
    Recombination on a dense (samples x loci) matrix.

    Each mating pair uses **one** ``get_breakpoints`` draw per pair. All siblings from
    that pair share the same crossover positions; within each segment, offspring rows
    partition inheritance between the two parents (rotating assignment along segments)
    so segment sources do not overlap redundantly across siblings the way independent
    per-offspring draws would.

    ``total_bp`` sums ``num_bp`` from that single draw per couple.
    """
    num_samples, genome_length = genotype_matrix.shape
    total_bp = 0

    num_pairs = num_samples // 2
    k = getattr(benchmark, "num_offspring_per_couple", 2)
    num_offspring = num_pairs * k

    offspring_matrix = np.zeros((num_offspring, genome_length), dtype=genotype_matrix.dtype)

    parent_indices = np.arange(num_samples)
    np.random.shuffle(parent_indices)

    offspring_idx = 0

    for i in range(0, num_samples - 1, 2):
        p1 = parent_indices[i]
        p2 = parent_indices[i + 1]

        if offspring_idx + k > num_offspring:
            return offspring_matrix, total_bp

        bp, num_bp = benchmark.get_breakpoints(bp_range, expected_crossovers)
        total_bp += num_bp

        _fill_shared_breakpoints_k_offspring(
            genotype_matrix,
            p1,
            p2,
            bp,
            k,
            offspring_idx,
            offspring_matrix,
        )
        offspring_idx += k

        if offspring_idx >= num_offspring:
            return offspring_matrix, total_bp

    return offspring_matrix, total_bp