covered(seg, I):      # seg fully covered by I
    l,r = seg.l, seg.r
    i=lower_bound(I, l)  # first interval with end>l
    cur=l
    while i<len(I) and cur<r:
        L,R = I[i].l, I[i].r
        if L>cur: return False
        cur = max(cur, R); i+=1
    return cur>=r

# Fast union of parent intervals for a node
parent_union(node):
    U=[]
    for (p, P) in node.parents: U += P
    return union_sets(U)

# Filter node.mutations to [l,r) in one pass (mutations sorted)
mut_slice(muts, l, r):
    # binary search to start, then linear until >=r
    i = lower_bound_pos(muts, l); out=[]
    while i<len(muts) and muts[i] < r: out.append(muts[i]); i+=1
    return out


# Attach offspring over seg=[L,R) using the shortest path:
# priority: reuse parents directly if possible; else structural add; else copy.
recombine_attach(node_id, offspring_id, L, R):
    seg=[L,R)
    stack=[node_id]
    visited=set()

    while stack:
        x = stack.pop()
        if x in visited: continue
        visited.add(x)

        nx = Nodes[x]
        Ux = parent_union(nx)             # O(k log k + k) small k expected
        Ix = intersect_sets(Ux, [seg])    # O(k)
        if not Ix:                        # no coverage here; push parents
            for (p, _) in nx.parents: stack.append(p)
            continue

        muts = mut_slice(nx.mutations, L, R)  # O(m_seg + log m)

        # Fast path: fully covered by parents and no local muts -> wire parents->offspring
        if not muts and covered(seg, Ux):
            for (p, P) in nx.parents:
                parts = intersect_sets(P, [seg])
                for part in parts: Edges.append((p, offspring_id, part))
            # push parents once to ensure upstream exposure
            for (p, _) in nx.parents: stack.append(p)
            continue

        # Mixed or has muts:
        if ADDING_NODES:
            # Try to cover seg directly from parents; if not, fall back to via self overlaps
            new_idx = alloc_id()
            new_node = Node(new_idx)
            new_node.mutations = muts[:]        
            added = False
            for (p, P) in nx.parents:
                if covered(seg, P):
                    new_node.parents.append((p, [seg]))
                    Edges.append((p, new_idx, seg))
                    added = True
            if not added:
                for part in Ix:        
                    new_node.parents.append((x, [part]))
                    Edges.append((x, new_idx, part))
            Edges.append((new_idx, offspring_id, seg))

        elif COPYING_MUTATIONS:
            # Cheap path: copy muts and link any actual covering intervals directly
            for m in muts: add_mut(offspring_id, m)
            for (p, P) in nx.parents:
                for part in intersect_sets(P, [seg]):
                    Edges.append((p, offspring_id, part))