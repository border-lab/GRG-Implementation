import networkx as nx
import matplotlib.pyplot as plt
import statistics
def plot_grg(grg, figsize=(14, 10)):
    G = nx.DiGraph()
    for node in grg.nodes:
        muts = [f"m{m}" for m in node.mut if m is not None]
        label = f"{node.n}\n" + ",".join(muts) if muts else f"{node.n}"
        G.add_node(node, label=label)

    for node in grg.nodes:
        for child in node.children:
            G.add_edge(node, child)

    endpoint_nodes = set(grg.haplotype_endpoints.values())

    height = {n: 0 for n in G.nodes}
    
    try:
        topo_order = list(nx.topological_sort(G))
    except nx.NetworkXUnfeasible:
        topo_order = list(G.nodes)

    for n in reversed(topo_order):
        if n in endpoint_nodes:
            height[n] = 0
        else:
            children = list(G.successors(n))
            if children:
                height[n] = max(height[c] for c in children) + 1
            else:
                height[n] = 1 

    layers = {}
    for n, h in height.items():
        layers.setdefault(h, []).append(n)
    
    max_height = max(layers.keys()) if layers else 0
    
    pos_x = {}
    min_separation = 4.0
    
    for h, nodes in layers.items():
        width = (len(nodes) - 1) * min_separation
        start_x = -width / 2
        for i, n in enumerate(nodes):
            pos_x[n] = start_x + i * min_separation
    
    iterations = 20
    for _ in range(iterations):
        directions = [range(max_height + 1), range(max_height, -1, -1)]
        
        for h_range in directions:
            for h in h_range:
                if h not in layers: continue
                nodes = layers[h]
            
                for n in nodes:
                    neighbors = list(G.successors(n)) + list(G.predecessors(n))
                    valid_neighbors = [nbr for nbr in neighbors if nbr in pos_x]
                    if valid_neighbors:
                        avg_x = statistics.mean([pos_x[nbr] for nbr in valid_neighbors])
                        pos_x[n] = avg_x
                nodes.sort(key=lambda n: pos_x[n])
                for i in range(len(nodes) - 1):
                    u, v = nodes[i], nodes[i+1]
                    if pos_x[v] - pos_x[u] < min_separation:
                        mid = (pos_x[u] + pos_x[v]) / 2
                        pos_x[u] = mid - min_separation / 2
                        pos_x[v] = mid + min_separation / 2

                if nodes:
                    current_center = (pos_x[nodes[0]] + pos_x[nodes[-1]]) / 2
                    for n in nodes:
                        pos_x[n] -= current_center

    pos = {}
    for n, h in height.items():
        pos[n] = (pos_x[n], h * 5.0) 

    plt.figure(figsize=figsize)
    
    node_colors = ["lightcoral" if n in endpoint_nodes else "skyblue" for n in G.nodes]
    node_edges = ["red" if n in endpoint_nodes else "black" for n in G.nodes]

    nx.draw(
        G, pos,
        with_labels=False,
        node_size=2200,
        node_color=node_colors,
        edgecolors=node_edges,
        edge_color="gray",
        arrowsize=15,
        width=1.5,
        alpha=0.9
    )

    labels = {n: data["label"] for n, data in G.nodes(data=True)}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=9, font_weight="bold")

    plt.title("Genome Representation Graph (GRG)", fontsize=14, fontweight="bold")
    plt.axis("off")
    plt.tight_layout()
    plt.show()
