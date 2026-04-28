"""
grg_recombination.py

Optimized core module containing the Non-duplication GRG recombination algorithm.
"""

import numpy as np
import bisect

# def get_breakpoints(N, expected_crossovers=1.5):
#     num_bp = np.random.poisson(expected_crossovers)
#     if num_bp == 0:
#         return np.array([], dtype=int)
#     num_bp = min(num_bp, N - 1)
#     bp = np.random.choice(range(1, N), size=num_bp, replace=False)
#     bp.sort()
#     return bp

# def recombination_intervals(h1, h2, N, expected_crossovers=1.5):
#     bp = get_breakpoints(N, expected_crossovers)
#     start = np.random.binomial(1, 0.5, 1)[0]
#     parents = [h1, h2]
    
#     segments = []
#     for i, K in enumerate(bp):
#         segments.append((parents[(start + i) % 2], K))
#     segments.append((parents[(start + len(bp)) % 2], N))
#     return segments

def get_breakpoints(bp_range, expected_crossovers=1.5):
    num_bp = np.random.poisson(expected_crossovers)
    if num_bp == 0:
        return np.array([], dtype=int)

    low, high = bp_range            # inclusive or exclusive? assuming inclusive high here
    length = high - low + 1
    num_bp = min(num_bp, length)

    # Oversample a bit, then deduplicate
    k = num_bp * 3                  # oversampling factor; small, since num_bp is tiny
    candidates = np.random.randint(low, high + 1, size=k)
    unique = np.unique(candidates)

    if unique.size < num_bp:
        # Very unlikely with small num_bp, but handle just in case
        extra_needed = num_bp - unique.size
        extra = np.random.randint(low, high + 1, size=extra_needed * 3)
        unique = np.unique(np.concatenate([unique, extra]))

    bp = np.sort(unique[:num_bp])
    return bp

def recombination_intervals(h1, h2, bp, N):
    """
    Returns list of segments: [(source_parent_id, end_coord), ...]
    """
    # bp = self.get_breakpoints(N)
    # start = np.random.binomial(1, 0.5, 1)[0]
    start = h1
    parents = [h1, h2]
    
    segments = []
    for i, K in enumerate(bp):
        segments.append((parents[(start + i) % 2], K))
    segments.append((parents[(start + len(bp)) % 2], N))
    return segments

class NonDuplicationRecombination:
    """
    Non-duplication GRG recombination algorithm.
    
    Key operations:
    - Path Compression: Skip nodes with no relevant mutations
    - Pruning: Stop when ancestry is disjoint from query interval
    - Bubble Insertion: Split nodes when partial mutation inheritance is needed
    """

    debug_mode = False  # Set to True to enable debug prints
    
    # def __init__(self, grg):
    #     self.grg = grg
    #     self.genome_length = grg.bp_range[1]  # l_max
    #     self.original_bp_range = grg.bp_range  # Store original bp_range for reference
    #     self._span_cache = {}               # NEW: Caches the full ancestor span
    #     self._ancestral_coverage_cache = {} # NEW: Caches the ancestral coverage
    #     self._mutation_cache = {}
    #     self.NEGATIVE_NODE_IDS = []
    #     self._modified_nodes = set()        # Track which nodes were modified for selective cache invalidation
    #     self._pending_bubbles = []          # Queue of bubble operations to defer until after traversal
    
    def __init__(self, grg):
        self.grg = grg
        self.genome_length = grg.bp_range[1]
        self.original_bp_range = grg.bp_range
        self.NEGATIVE_NODE_IDS = []
        self._modified_nodes = set()
        self._pending_bubbles = []
        self._pending_sample_removals = set()

        # JIT Caches: These start empty and only populate for visited nodes
        self._mutation_cache = {}  # Stores [(id, pos), ...]
        self._pos_cache = {}       # Stores [pos, ...] for bisect
        
        # Keep these as lists/arrays for O(1) DP access
        self.span_cache = [False] * self.grg.num_nodes
        self.anc_cov_cache = [False] * self.grg.num_nodes

    # def _prefetch_all_mutations(self):
    #     self._mutation_cache = {i: [] for i in range(self.grg.num_nodes)}
    #     self._pos_cache = {} # NEW: Raw float lists for bisect

    #     for node_id, mut_id in self.grg.get_node_mutation_pairs():
    #         mut = self.grg.get_mutation_by_id(mut_id)
    #         self._mutation_cache[node_id].append((mut_id, mut.position))
            
    #     for node_id in self._mutation_cache:
    #         if self._mutation_cache[node_id]:
    #             self._mutation_cache[node_id].sort(key=lambda x: x[1])
    #             # Lock in a pure list of floats for O(log M) searching
    #             self._pos_cache[node_id] = [m[1] for m in self._mutation_cache[node_id]]

    # def _get_node_mutations(self, node_id):
    #     """Get mutation positions for a node (cached, sorted by position)."""
    #     if node_id not in self._mutation_cache:
    #         mut_ids = self.grg.get_mutations_for_node(node_id)
    #         mutations = []
    #         for mut_id in mut_ids:
    #             mut = self.grg.get_mutation_by_id(mut_id)
    #             mutations.append((mut_id, mut.position))
    #         # Sort by position for fast interval queries using binary search
    #         mutations.sort(key=lambda x: x[1])
    #         self._mutation_cache[node_id] = mutations
    #     return self._mutation_cache[node_id]

    def _get_node_mutations(self, node_id):
        """Get mutation positions for a node (Loaded JIT and cached)."""
        if node_id not in self._mutation_cache:
            # Only fetch from C++ if not in cache
            mut_ids = self.grg.get_mutations_for_node(node_id)
            mutations = []
            positions = []
            
            for mut_id in mut_ids:
                mut = self.grg.get_mutation_by_id(mut_id)
                mutations.append((mut_id, mut.position))
                positions.append(mut.position)
            
            # Sort once and cache both the full data and the raw positions
            # Use the position list specifically for fast bisecting
            combined = sorted(mutations, key=lambda x: x[1])
            self._mutation_cache[node_id] = combined
            self._pos_cache[node_id] = [m[1] for m in combined]
            
        return self._mutation_cache[node_id]
    
    def _get_mutations_in_interval(self, node_id, L, R):
        """Get mutations on node_id that fall within [L, R) using binary search (O(log N) time)."""
        node_muts = self._get_node_mutations(node_id)
        if not node_muts:
            return []
        
        # Extract positions for binary search
        positions = [pos for _, pos in node_muts]
        
        # Find the range [L, R) using binary search
        left_idx = bisect.bisect_left(positions, L)
        right_idx = bisect.bisect_left(positions, R)
        
        return node_muts[left_idx:right_idx]
    
    def _has_connected_descendant(self, node_id, connected, visited=None):
        """True if any descendant of node_id (via down edges) is in connected."""
        if visited is None:
            visited = set()
        if node_id in visited:
            return False
        visited.add(node_id)
        for child in self.grg.get_down_edges(node_id):
            if child in connected:
                return True
            if self._has_connected_descendant(child, connected, visited):
                return True
        return False

    def _has_connected_ancestor(self, node_id, connected, visited=None):
        """True if any ancestor of node_id (via up edges) is in connected."""
        if visited is None:
            visited = set()
        if node_id in visited:
            return False
        visited.add(node_id)
        for parent in self.grg.get_up_edges(node_id):
            if parent in connected:
                return True
            if self._has_connected_ancestor(parent, connected, visited):
                return True
        return False

    # def _get_ancestral_coverage(self, node_id):
    #     """Compute the ancestral interval coverage Iu (mutations in ancestors ONLY)."""
    #     # 1. Check cache
    #     if node_id in self._ancestral_coverage_cache:
    #         return self._ancestral_coverage_cache[node_id]

    #     parents = self.grg.get_up_edges(node_id)
    #     if not parents:
    #         self._ancestral_coverage_cache[node_id] = None
    #         return None
            
    #     min_pos, max_pos = float('inf'), float('-inf')
    #     visited = set() 
        
    #     for parent in parents:
    #         parent_span = self._get_node_and_ancestor_span(parent, visited)
    #         if parent_span:
    #             min_pos = min(min_pos, parent_span[0])
    #             max_pos = max(max_pos, parent_span[1])
        
    #     # 2. Format result and cache it
    #     result = None if min_pos == float('inf') else (min_pos, max_pos + 1)
        
    #     self._ancestral_coverage_cache[node_id] = result
    #     return result

    # def _get_node_and_ancestor_span(self, node_id, visited):
    #     """Recursive helper: gets span of mutations for this node and its ancestors."""
    #     # 1. Check cache
    #     if node_id in self._span_cache:
    #         return self._span_cache[node_id]
            
    #     if node_id in visited:
    #         return None
    #     visited.add(node_id)
        
    #     min_pos, max_pos = float('inf'), float('-inf')
        
    #     # 2. Check this node's mutations
    #     node_muts = self._get_node_mutations(node_id)
    #     if node_muts:
    #         min_pos = min(min_pos, min(pos for _, pos in node_muts))
    #         max_pos = max(max_pos, max(pos for _, pos in node_muts))
        
    #     # 3. Recurse to parents
    #     for parent in self.grg.get_up_edges(node_id):
    #         anc_span = self._get_node_and_ancestor_span(parent, visited)
    #         if anc_span:
    #             min_pos = min(min_pos, anc_span[0])
    #             max_pos = max(max_pos, anc_span[1])
                
    #     # 4. Format result and cache it
    #     result = None if min_pos == float('inf') else (min_pos, max_pos)
        
    #     self._span_cache[node_id] = result
    #     return result

    def _get_mutation_range(self, node_id, L, R):
        """Returns (start_idx, end_idx) for mutations in [L, R) in O(log M)."""
        # Ensure data is loaded JIT
        self._get_node_mutations(node_id)
            
        positions = self._pos_cache[node_id]
        if not positions:
            return 0, 0
            
        return bisect.bisect_left(positions, L), bisect.bisect_left(positions, R)

    def _get_node_and_ancestor_span(self, node_id):
        """Lazy recursive helper using array memoization."""
        # 1. Instant Array Cache Hit
        if self.span_cache[node_id] is not False:
            return self.span_cache[node_id]
            
        min_pos, max_pos = float('inf'), float('-inf')
        
        # 2. Check this node's mutations
        node_muts = self._get_node_mutations(node_id)
        if node_muts:
            min_pos = node_muts[0][1]
            max_pos = node_muts[-1][1]
        
        # 3. Recurse to parents (only runs once per node)
        for parent in self.grg.get_up_edges(node_id):
            anc_span = self._get_node_and_ancestor_span(parent)
            if anc_span:
                min_pos = min(min_pos, anc_span[0])
                max_pos = max(max_pos, anc_span[1])
                
        # 4. Format result and lock it in the array
        result = None if min_pos == float('inf') else (min_pos, max_pos)
        self.span_cache[node_id] = result
        return result

    def _get_ancestral_coverage(self, node_id):
        """Compute Iu lazily using array memoization."""
        if self.anc_cov_cache[node_id] is not False:
            return self.anc_cov_cache[node_id]

        parents = self.grg.get_up_edges(node_id)
        if not parents:
            self.anc_cov_cache[node_id] = None
            return None
            
        min_pos, max_pos = float('inf'), float('-inf')
        
        for parent in parents:
            p_span = self._get_node_and_ancestor_span(parent)
            if p_span:
                min_pos = min(min_pos, p_span[0])
                max_pos = max(max_pos, p_span[1])
                
        result = None if min_pos == float('inf') else (min_pos, max_pos + 1)
        self.anc_cov_cache[node_id] = result
        return result   
    
    def _extract_bubble(self, node_id, relevant_mut_ids, offspring_id, interval):
        """
        Queue a bubble node creation to be deferred until after traversal.
        
        Instead of immediately modifying the graph, this method queues the bubble
        operation to be applied after traversal completes. This prevents mid-traversal
        cache invalidations and graph disruptions.
        
        Creates a new node v that:
        - Contains the relevant mutations (moved from node_id)
        - Becomes a parent of node_id
        - Becomes a parent of the offspring
        
        Args:
            node_id: The node to split
            relevant_mut_ids: Mutation IDs to move to the bubble
            offspring_id: The offspring node to connect
            interval: The interval [L, R) being inherited
            
        Returns:
            The bubble node ID (created immediately for graph structure)
        """

        # Create new bubble node immediately (needed for edge connections)
        bubble_id = self.grg.make_node()

        # Expand Lazy arrays for the new node safely
        while len(self.span_cache) <= bubble_id:
            self.span_cache.append(False)
            self.anc_cov_cache.append(False)         
        
        # Add connections immediately (needed for graph structure)
        self.grg.connect(bubble_id, node_id)
        self.grg.connect(bubble_id, -offspring_id)
        
        # Queue mutations to be moved later (deferred)
        # This avoids cascading cache invalidations during traversal
        self._pending_bubbles.append({
            'node_id': node_id,
            'bubble_id': bubble_id,
            'relevant_mut_ids': relevant_mut_ids
        })
        
        # Track which nodes were modified for selective cache invalidation (to be done after traversal)
        self._modified_nodes.add(node_id)        # Will lose mutations
        self._modified_nodes.add(bubble_id)      # New bubble node
        # Also track parents since node_id's composition will change
        for parent in self.grg.get_up_edges(node_id):
            self._modified_nodes.add(parent)
        
        # Don't invalidate caches yet - let traversal continue on stable graph
        # They'll be cleared after all bubbles are applied
        
        return bubble_id
    
    # def _recurse_attach(self, node_id, offspring_id, L, R, visited=None, connected=None):
    #     """
    #     Recursively attach ancestry to offspring for interval [L, R).
        
    #     Algorithm 2: RecurseAttach
        
    #     Handles:
    #     - Full Coverage (I ⊆ Iu): Ancestors fully cover query
    #     - Disjoint (I ∩ Iu = ∅): Pruning - stop traversal
    #     - Partial Overlap: Decomposition - recurse on intersection
        
    #     Args:
    #         node_id: Current ancestor node
    #         offspring_id: The new offspring node
    #         L, R: Query interval [L, R)
    #         visited: Set of already visited nodes to avoid cycles
    #         connected: Set of nodes already connected to offspring (avoid duplicates)
    #     """
    #     if L >= R:
    #         return
        
    #     if visited is None:
    #         visited = set()
    #     if connected is None:
    #         connected = set()
        
    #     if node_id in visited:
    #         return
    #     visited.add(node_id)
        
    #     # # Get mutations and coverage info
    #     # all_muts = self._get_node_mutations(node_id)
    #     # relevant_muts = self._get_mutations_in_interval(node_id, L, R)
        
    #     # all_mut_ids = set(mut_id for mut_id, _ in all_muts)
    #     # relevant_mut_ids = set(mut_id for mut_id, _ in relevant_muts)

    #     # # Check mutation relevance
    #     # has_all_relevant = (relevant_mut_ids == all_mut_ids) and len(all_mut_ids) > 0
    #     # has_no_relevant = len(relevant_mut_ids) == 0
    #     # has_partial_relevant = len(relevant_mut_ids) > 0 and relevant_mut_ids != all_mut_ids

    #     # Get mutations and coverage info
    #     all_muts = self._get_node_mutations(node_id)
    #     relevant_muts = self._get_mutations_in_interval(node_id, L, R)
        
    #     num_all = len(all_muts)
    #     num_rel = len(relevant_muts)

    #     # NEW: FAST O(1) mathematical checks instead of O(M) set creations
    #     has_all_relevant = (num_all > 0) and (num_rel == num_all)
    #     has_no_relevant = (num_rel == 0)
    #     has_partial_relevant = (num_rel > 0) and (num_rel < num_all)
        
    #     # Get ancestral coverage
    #     Iu = self._get_ancestral_coverage(node_id)
        
    #     if self.debug_mode:
    #         print(f"Visiting node {node_id}: all_muts={all_muts}, relevant_muts={relevant_muts}, Iu={Iu}")
        
    #     # Determine coverage status
    #     if Iu is None:
    #         if self.debug_mode:
    #             print("No ancestral mutations - treating as root node")
    #         if has_all_relevant:
    #             if self.debug_mode:
    #                 print("Root Node Case 1: All relevant mutations - connect directly")
    #             # if (node_id not in connected and
    #             #     not self._has_connected_descendant(node_id, connected)):
    #             if node_id not in connected:
    #                 self.grg.connect(node_id, -offspring_id)
    #                 connected.add(node_id)

    #                 # Queue sample removal instead of crossing into C++ immediately
    #                 self._pending_sample_removals.add(node_id)

    #                 # Remove sample status if connecting an original sample node to avoid confusion with new offspring samples
    #                 # current_samples = list(self.grg.get_sample_nodes())
    #                 # if node_id in current_samples:
    #                 #     current_samples.remove(node_id)
    #                 #     self.grg.set_samples(current_samples)

    #             elif self.debug_mode:
    #                 print(f"Node {node_id} already connected or has connected relatives, skipping direct connection")
            
    #         if has_partial_relevant:
    #             if self.debug_mode:
    #                 print("Root Node Case 2: Partial relevant mutations - create bubble and connect")
                
    #             # Generate the list of IDs ONLY if a bubble is actually needed
    #             rel_mut_ids = [mut_id for mut_id, _ in relevant_muts]
    #             bubble_id = self._extract_bubble(node_id, rel_mut_ids, offspring_id, (L, R))
    #             connected.add(bubble_id)
                
    #             # bubble_id = self._extract_bubble(node_id, list(relevant_mut_ids), 
    #             #                                         offspring_id, (L, R))
    #             # connected.add(bubble_id)
    #             # connected.add(node_id)
            
    #         return # Stop - lineage fully resolved

    def _recurse_attach(self, node_id, offspring_id, L, R, visited=None, connected=None):
        if L >= R: return
        
        if visited is None: visited = set()
        if connected is None: connected = set()
        
        if node_id in visited: return
        visited.add(node_id)
        
        # 1. Get ONLY the range indices (O(log M))
        left, right = self._get_mutation_range(node_id, L, R)
        num_rel = right - left
        num_all = len(self._pos_cache.get(node_id, []))

        has_all_relevant = (num_all > 0) and (num_rel == num_all)
        has_no_relevant = (num_rel == 0)
        has_partial_relevant = (num_rel > 0) and (num_rel < num_all)
        
        Iu = self._get_ancestral_coverage(node_id)
        
        # Scenario: Root / Terminal Node
        if Iu is None:
            if has_all_relevant:
                if node_id not in connected:
                    self.grg.connect(node_id, -offspring_id)
                    connected.add(node_id)
                    self._pending_sample_removals.add(node_id)
            
            if has_partial_relevant:
                # ONLY NOW do we extract the IDs for the bubble
                rel_mut_ids = [m[0] for m in self._mutation_cache[node_id][left:right]]
                bubble_id = self._extract_bubble(node_id, rel_mut_ids, offspring_id, (L, R))
                connected.add(bubble_id)
            return
        
        # Check if disjoint (Pruning - Scenario 2)
        # I ∩ Iu = ∅
        ancestral_disjoint = R <= Iu[0] or L >= Iu[1]

        # Check interval coverage
        full_coverage = (Iu[0] >= L and Iu[1] <= R)  # Iu ⊆ I (node is fully consumed)
        
        if has_all_relevant and full_coverage:
            if self.debug_mode:
                print("Case 1: Full coverage with all relevant mutations - connect directly")
            # if (node_id not in connected and
            #     not self._has_connected_descendant(node_id, connected)):
            if node_id not in connected:
                self.grg.connect(node_id, -offspring_id)
                connected.add(node_id)

                # Queue sample removal instead of crossing into C++ immediately
                self._pending_sample_removals.add(node_id)

                # Remove sample status if connecting an original sample node to avoid confusion with new offspring samples
                # current_samples = list(self.grg.get_sample_nodes())
                # if node_id in current_samples:
                #     current_samples.remove(node_id)
                #     self.grg.set_samples(current_samples)

            elif self.debug_mode:
                print(f"Node {node_id} already connected or has connected relatives, skipping direct connection")
            return  # Stop - lineage fully resolved
        
        if has_no_relevant and ancestral_disjoint:
            if self.debug_mode:
                print("Case 2: No relevant mutations and disjoint ancestry - prune")
            # End search no relevant mutations to be found
            return
        
        if has_no_relevant and not ancestral_disjoint:
            if self.debug_mode:
                print("Case 3: No relevant mutations but not disjoint ancestry - path compression")
            # Path Compression - bypass this node
            # Recurse on parents
            parents = self.grg.get_up_edges(node_id)
            
            # If there's partial overlap with ancestry, we need to adjust the interval for the parents to the intersection of [L, R) and Iu
            newL = max(L, Iu[0]) 
            newR = min(R, Iu[1]) 
            L, R = newL, newR
            if L >= R: #avoids useless parent calls on degenerate overlap at boundaries 
                    return

            for parent in parents:
                self._recurse_attach(parent, offspring_id, L, R, visited, connected)
            return

        if has_partial_relevant or has_all_relevant:
            if self.debug_mode:
                print("Case 4: Partial/full relevant mutations - need to create bubble and maybe recurse upwards")
            # Cases where we have relevant mutations (partial or full) but not full coverage - need to create bubble and maybe recurse upwards
            
            # Generate the list of IDs ONLY if a bubble is actually needed
            rel_mut_ids = [m[0] for m in self._mutation_cache[node_id][left:right]]
            bubble_id = self._extract_bubble(node_id, rel_mut_ids, offspring_id, (L, R))
            connected.add(bubble_id)
            
            # bubble_id = self._extract_bubble(node_id, list(relevant_mut_ids), 
            #                                         offspring_id, (L, R))
            # connected.add(bubble_id)
            # connected.add(node_id)  # Treat split node as covered to avoid connecting ancestors

            if not ancestral_disjoint:
                # If there's partial overlap with ancestry, we need to adjust the interval for the parents to the intersection of [L, R) and Iu
                newL = max(L, Iu[0])
                newR = min(R, Iu[1])
                L, R = newL, newR
                if L >= R: #avoids useless parent calls on degenerate overlap at boundaries 
                    return

                parents = self.grg.get_up_edges(node_id)
                for parent in parents:
                    self._recurse_attach(parent, offspring_id, L, R, visited, connected)
            
            return
    
    # def _apply_pending_bubbles(self):
    #     """Apply all deferred bubble modifications after traversal completes.
        
    #     This batches all mutation movements and ensures the graph structure
    #     remains stable during traversal, avoiding cascading cache invalidations.
    #     """
    #     for bubble_op in self._pending_bubbles:
    #         node_id = bubble_op['node_id']
    #         bubble_id = bubble_op['bubble_id']
    #         relevant_mut_ids = bubble_op['relevant_mut_ids']
            
    #         # Move relevant mutations from node to bubble (Algorithm 3: v.M <- Mrel; u.M <- u.M \ Mrel)
    #         for mut_id in relevant_mut_ids:
    #             mut = self.grg.get_mutation_by_id(mut_id)
                
    #             # Add mutation to bubble node
    #             self.grg.add_mutation(mut, bubble_id)
    #             # Remove mutation from original node
    #             self.grg.remove_mutation(mut_id, node_id)
        
    #     # Clear the queue
    #     self._pending_bubbles.clear()

    def _apply_pending_bubbles(self):
        """Apply all deferred bubble modifications and sample removals after traversal."""
        for bubble_op in self._pending_bubbles:
            node_id = bubble_op['node_id']
            bubble_id = bubble_op['bubble_id']
            relevant_mut_ids = bubble_op['relevant_mut_ids']
            
            for mut_id in relevant_mut_ids:
                mut = self.grg.get_mutation_by_id(mut_id)
                self.grg.add_mutation(mut, bubble_id)
                self.grg.remove_mutation(mut_id, node_id)
        
        self._pending_bubbles.clear()

        # NEW: Process pending sample removals in one single fast batch
        if self._pending_sample_removals:
            current_samples = set(self.grg.get_sample_nodes())
            current_samples.difference_update(self._pending_sample_removals)
            self.grg.set_samples(list(current_samples))
            self._pending_sample_removals.clear()
    
    # def _clear_modified_caches(self):
    #     """Clear only the caches for nodes that were modified during recombination.
        
    #     This is much more efficient than clearing all caches, as most nodes in a large
    #     graph are unmodified and their cached data remains valid for the next offspring.
    #     """
    #     for node_id in self._modified_nodes:
    #         self._mutation_cache.pop(node_id, None)
    #         self._span_cache.pop(node_id, None)
    #         self._ancestral_coverage_cache.pop(node_id, None)
    #     self._modified_nodes.clear()
    #     self._precompute_spans_iteratively()  # Recompute spans for modified nodes and their ancestors

    def _clear_modified_caches(self):
        for node_id in self._modified_nodes:
            self._mutation_cache.pop(node_id, None)
            self._pos_cache.pop(node_id, None) # Clear the position cache too!
            if node_id < len(self.span_cache):
                self.span_cache[node_id] = False
                self.anc_cov_cache[node_id] = False
        self._modified_nodes.clear()
    
    def recombine(self, haplotype_A, haplotype_B, breakpoint):
        """
        Generate offspring through recombination.
        
        Creates a new offspring node inheriting:
        - [0, breakpoint) from haplotype_A
        - [breakpoint, genome_length) from haplotype_B
        
        Args:
            haplotype_A: Node ID of first parent
            haplotype_B: Node ID of second parent  
            breakpoint: Crossover position
            
        Returns:
            Node ID of the new offspring
        """
        # Ensure pending bubbles are cleared from previous operations
        self._pending_bubbles.clear()
        
        # Create offspring node
        offspring_id = self.grg.make_node(negative=True)
        
        # Track connected nodes to avoid duplicate edges
        connected = set()
        
        # Inherit from haplotype_A for [0, breakpoint)
        self._recurse_attach(haplotype_A, offspring_id, 0, breakpoint, 
                            visited=None, connected=connected)
        
        # Inherit from haplotype_B for [breakpoint, genome_length)
        self._recurse_attach(haplotype_B, offspring_id, breakpoint, self.genome_length,
                            visited=None, connected=connected)
        
        # Apply all deferred bubble modifications after traversal completes
        self._apply_pending_bubbles()

        # NEW: Clear caches and rebuild DP arrays for the next generation
        self._clear_modified_caches()
        
        #self.grg.bp_range = self.original_bp_range  # Ensure bp_range is updated to reflect the full genome length
        
        if offspring_id not in self.NEGATIVE_NODE_IDS:
            self.NEGATIVE_NODE_IDS.append(offspring_id)
        return -(self.NEGATIVE_NODE_IDS.index(offspring_id) + 1)
    
    def recombine_multi(self, segments):
        """
        Generate offspring from multiple segments.
        
        Args:
            segments: List of (parent_node_id, interval_end) tuples
                     where intervals are implicit from previous end
                     
        Returns:
            Node ID of the new offspring
        """
        # Ensure pending bubbles are cleared from previous operations
        self._pending_bubbles.clear()
        
        # Create offspring node
        offspring_id = self.grg.make_node(negative=True)
        
        # Track connected nodes to avoid duplicate edges (maintained across all segments)
        connected = set()
        
        # Process each segment
        start = 0
        for parent_id, end in segments:
            if end > start:
                if self.debug_mode:
                    print("BREAK")
                # Reuse connected set across segments to avoid duplicate connections
                self._recurse_attach(parent_id, offspring_id, start, end,
                                    visited=None, connected=connected)
            start = end
        
        # Apply all deferred bubble modifications after traversal completes
        # This prevents mid-traversal cache invalidations and keeps the graph stable
        self._apply_pending_bubbles()

        # NEW: Clear caches and rebuild DP arrays for the next generation
        self._clear_modified_caches()
        
        #self.grg.bp_range = self.original_bp_range  # Ensure bp_range is updated to reflect the full genome length
        
        if offspring_id not in self.NEGATIVE_NODE_IDS:
            self.NEGATIVE_NODE_IDS.append(offspring_id)
        return -(self.NEGATIVE_NODE_IDS.index(offspring_id) + 1)

def simulate_grg_recombination(recomb, bp_range, N):
    # recombiner = NonDuplicationRecombination(grg)
    
    samples = np.array(recomb.grg.get_sample_nodes())
    np.random.shuffle(samples)

    new_offspring_ids = []
    
    for i in range(0, len(samples) - 1, 2):
        p1 = samples[i]
        p2 = samples[i+1]
        
        for j in range(2):
            bp = get_breakpoints(bp_range, expected_crossovers=1.5)
            segments = recombination_intervals(p1, p2, bp, N)
            
            # Clear only caches for nodes that were modified in the previous offspring.
            # This preserves cache data for the ~95% of nodes that weren't touched,
            # dramatically reducing recomputation compared to clearing all caches.
            recomb._clear_modified_caches()
            
            offspring_id = recomb.recombine_multi(segments)
            raw_id = recomb.NEGATIVE_NODE_IDS[abs(offspring_id) - 1]
            new_offspring_ids.append(raw_id)
            
    new_offspring_ids.sort()
    recomb.grg.set_samples(new_offspring_ids)
        
    return new_offspring_ids