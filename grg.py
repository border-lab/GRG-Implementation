from collections import defaultdict

import numpy as np

class Node:
    def __init__(self, n: int, mut):
        self.n = n
        self.mut = mut if isinstance(mut, tuple) else (mut,)
        self.parents = {}
        self.children = {}

    def __repr__(self):
        return f'Node("{self.n}, {self.mut}")'

    def _normalize(self, intervals):
        """Sort and merge overlapping/touching [l,r) intervals."""
        if not intervals:
            return []
        intervals = sorted(intervals, key=lambda x: x[0])
        out = [intervals[0][:]]
        for l, r in intervals[1:]:
            L, R = out[-1]
            if l <= R:
                out[-1][1] = max(R, r)
            else:
                out.append([l, r])
        return out

    def connect(self, child: "Node", interval):
        if not (isinstance(interval, (list, tuple)) and len(interval) == 2):
            raise ValueError("interval must be [l, r]")

        l, r = interval
        if not (isinstance(l, (int, float)) and isinstance(r, (int, float)) and l < r):
            raise ValueError("interval must be numeric [l, r) with l < r")

        # Update children map
        lst = self.children.get(child, [])
        lst.append([l, r])
        self.children[child] = self._normalize(lst)

        # Update child's parents map
        plst = child.parents.get(self, [])
        plst.append([l, r])
        child.parents[self] = child._normalize(plst)

    def verify(self):
        """Type / structural checks for this node."""
        if not isinstance(self.n, int):
            raise TypeError("Node ID 'n' must be an integer.")
        if not isinstance(self.mut, tuple):
            raise TypeError("Node 'mut' attribute must be a tuple.")

        if not isinstance(self.parents, dict):
            raise TypeError("Node 'parents' must be a dict.")
        if not isinstance(self.children, dict):
            raise TypeError("Node 'children' must be a dict.")

        for p, ivals in self.parents.items():
            if not isinstance(p, Node):
                raise TypeError("Parents must be Node instances.")
            if not isinstance(ivals, list):
                raise TypeError("Parent intervals must be a list.")
        for c, ivals in self.children.items():
            if not isinstance(c, Node):
                raise TypeError("Children must be Node instances.")
            if not isinstance(ivals, list):
                raise TypeError("Child intervals must be a list.")
        return True


class GRG:
    def __init__(
        self,
        nodes: set[Node],
        haplotype_endpoints,          # dict[hap_id] -> terminal Node
        haplotype_to_individual_map,  # dict[hap_id] -> individual_id
        samples: list[str] = None,
    ):
        self.nodes = nodes
        self.mutations = [
            m for node in nodes for m in node.mut if m is not None
        ]
        self.samples = samples if samples is not None else []
        self.haplotype_endpoints = haplotype_endpoints
        self.haplotype_to_individual_map = haplotype_to_individual_map
        self.verify()

    def verify(self):
        if not isinstance(self.nodes, set):
            raise TypeError("GRG 'nodes' attribute must be a set.")

        for node in self.nodes:
            node.verify()
            # Check children and ensure symmetric parent pointers
            for child, ivals in node.children.items():
                if not isinstance(child, Node):
                    raise TypeError(f"Node {node} has invalid child type {type(child)}.")
                if child not in self.nodes:
                    raise ValueError(f"Node {node} connects to child not in GRG.")
                # symmetry: this should appear in child's parents
                if node not in child.parents:
                    raise ValueError("Child.parents missing reciprocal link.")
            # Check parents symmetry
            for parent, ivals in node.parents.items():
                if not isinstance(parent, Node):
                    raise TypeError(f"Node {node} has invalid parent type {type(parent)}.")
                if parent not in self.nodes:
                    raise ValueError(f"Node {node} has parent not in GRG.")
                if node not in parent.children:
                    raise ValueError("Parent.children missing reciprocal link.")

        if self._has_cycles():
            raise ValueError("The graph contains a cycle and is not a valid DAG.")
        return True

    def _has_cycles(self):
        visiting = set()
        visited = set()

        def dfs(u: Node):
            visiting.add(u)
            visited.add(u)
            for v in u.children.keys():
                if v in visiting:
                    return True
                if v not in visited:
                    if dfs(v):
                        return True
            visiting.remove(u)
            return False

        for node in self.nodes:
            if node not in visited:
                if dfs(node):
                    return True
        return False

    def _get_all_parents(self):
        all_parents = {}
        for node in self.nodes:
            for parent in node.parents.keys():
                all_parents.setdefault(node, []).append(parent)

        all_children = set(all_parents.keys())
        roots = [node for node in self.nodes if node not in all_children]
        return all_parents, roots
    def get_individuals_to_mutations(self):
        all_parents, roots = self._get_all_parents()

        haplotype_mutations = {}
        for hap_id, terminal_node in self.haplotype_endpoints.items():
            muts = set()
            visited = set()
            stack = [terminal_node]
            while stack:
                current = stack.pop()
                if current in visited:
                    continue
                visited.add(current)
                muts.update(m for m in current.mut if m is not None)
                for parent in all_parents.get(current, []):
                    if parent not in visited:
                        stack.append(parent)
            haplotype_mutations[hap_id] = muts

        from collections import defaultdict
        individual_to_mut_sets = defaultdict(set)
        for hap_id, muts in haplotype_mutations.items():
            ind_id = self.haplotype_to_individual_map.get(hap_id)
            if ind_id:
                individual_to_mut_sets[ind_id].update(muts)

        return {ind: sorted(muts) for ind, muts in individual_to_mut_sets.items()}

    def get_mutation_to_individuals(self):
        individual_to_mutations = self.get_individuals_to_mutations()
        from collections import defaultdict
        mutation_to_ind_sets = defaultdict(set)
        for individual, mutations in individual_to_mutations.items():
            for mut in mutations:
                mutation_to_ind_sets[mut].add(individual)
        return {mut: sorted(list(inds)) for mut, inds in mutation_to_ind_sets.items()}

    # ---------- matrix representation ----------

    def verify_matrix(self, matrix, mutations_map):
        import numpy as np
        dims = matrix.shape
        if dims[0] != len(self.haplotype_endpoints):
            raise ValueError("Matrix doesn't have as many rows as haplotypes.")
        if dims[1] != len(mutations_map):
            raise ValueError("Matrix doesn't have as many columns as mutations.")
        if not np.all((matrix == 0) | (matrix == 1)):
            raise ValueError("Matrix breaches 0/1 rule.")
    def add_node(self, node: Node):
        self.nodes.add(node)

    def to_matrix(self):
        """
        Returns (matrix, mutation_order) where each row is a haplotype and
        each column is a mutation (1 if present on that haplotype).
        """
        import numpy as np

        all_mutations = sorted(set(m for node in self.nodes for m in node.mut if m is not None))
        mutation_to_col = {mut: idx for idx, mut in enumerate(all_mutations)}

        haplotype_ids = sorted(self.haplotype_endpoints.keys())
        matrix = np.zeros((len(haplotype_ids), len(all_mutations)), dtype=int)
        all_parents, _ = self._get_all_parents()

        for row_idx, hap_id in enumerate(haplotype_ids):
            terminal_node = self.haplotype_endpoints[hap_id]
            muts = set()
            visited = set()
            stack = [terminal_node]
            while stack:
                current = stack.pop()
                if current in visited:
                    continue
                visited.add(current)
                muts.update(m for m in current.mut if m is not None)
                for parent in all_parents.get(current, []):
                    stack.append(parent)
            for mut in muts:
                matrix[row_idx, mutation_to_col[mut]] = 1

        self.verify_matrix(matrix, all_mutations)
        return matrix, all_mutations

    # ---------- dot product over DAG ----------

    def dot(self, y):
        """
        Same semantics as before:
        sum_over_nodes( (sum y[m] over node.mut) * (#haplotype_endpoints_below_node) )
        now using children dict instead of adj set.
        """
        import numpy as np
        from collections import deque

        nodes_list = list(self.nodes)
        node_to_id = {node: i for i, node in enumerate(nodes_list)}
        n = len(nodes_list)

        node_mut_sum = np.zeros(n, dtype=np.float64)
        for i, node in enumerate(nodes_list):
            muts = [m for m in node.mut if m is not None]
            if muts:
                node_mut_sum[i] = np.sum(y[muts])

        # Build adjacency and indegree from children dict
        children = [[] for _ in range(n)]
        indegree = np.zeros(n, dtype=np.int64)
        for i, node in enumerate(nodes_list):
            for child in node.children.keys():
                j = node_to_id[child]
                children[i].append(j)
                indegree[j] += 1

        q = deque([i for i in range(n) if indegree[i] == 0])
        topo = []
        while q:
            u = q.popleft()
            topo.append(u)
            for v in children[u]:
                indegree[v] -= 1
                if indegree[v] == 0:
                    q.append(v)

        endpoint_ids = {node_to_id[node] for node in self.haplotype_endpoints.values()}
        endpoint_count = np.zeros(n, dtype=np.int64)
        for u in reversed(topo):
            count = 1 if u in endpoint_ids else 0
            for v in children[u]:
                count += endpoint_count[v]
            endpoint_count[u] = count

        total = np.sum(endpoint_count * node_mut_sum)
        return float(total)

