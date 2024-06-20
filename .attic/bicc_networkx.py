def biccs_nx(graph):
    def edge_stack_to_set(edge_stack):
        out_set = set()
        for es in edge_stack:
            for n in es:
                out_set.add(n)
        return out_set

    def next_child(stack_item):
        if not stack_item[3]:
            return None
        if stack_item[2] >= len(stack_item[3]):
            return None
        else:
            stack_item[2] += 1
            return stack_item[3][stack_item[2] - 1]

    # depth-first search algorithm to generate articulation points
    # and biconnected components
    visited = set()

    for n in graph.nodes:
        if n in visited:
            continue
        discovery = {n: 0}
        low = {n: 0}
        root_children = 0
        artic_points = []
        components = []
        visited.add(n)
        edge_stack = []
        # stack = [(start, start, iter(G[start]))]
        neighbors = graph[n].neighbors()
        stack = [[n, n, 0, neighbors]]
        while stack:
            parent = stack[-1][0]
            child = stack[-1][1]
            nn = next_child(stack[-1])

            # pdb.set_trace()
            if nn:
                if nn == parent:
                    continue
                if nn in visited:
                    if discovery[nn] <= discovery[child]:
                        edge_stack.append((child, nn))
                        low[child] = min(low[child], discovery[nn])
                else:
                    low[nn] = len(discovery)
                    discovery[nn] = low[nn]
                    visited.add(nn)
                    stack.append([child, nn, 0, graph[nn].neighbors()])
                    edge_stack.append((child, nn))
            elif nn is None:
                stack.pop()
                if len(stack) > 1:
                    if low[child] >= discovery[parent]:
                        artic_points.append(parent)
                        cut_point = edge_stack.index((parent, child))
                        comp = edge_stack_to_set(edge_stack[cut_point:])
                        components.append(comp)
                        edge_stack = edge_stack[:cut_point]
                        # components.append(edge_stack_to_set(edge_stack[edge_stack.index((parent, child)):]))
                    low[parent] = min(low[parent], low[child])
                elif stack:
                    root_children += 1
                    components.append(
                        edge_stack_to_set(edge_stack[edge_stack.index((parent, child)) :])
                    )

        if root_children > 1:
            artic_points.append[n]

    return components, artic_points

    def bi_cc_rec(
        self, n_id, parent, low, disc, stack, node_ids, all_bi_cc, artic_points, disc_time
    ):
        # Count of children in current node
        u = node_ids[n_id]
        children = 0
        # Initialize discovery time and low value
        disc[u] = disc_time[0]
        low[u] = disc_time[0]
        disc_time[0] += 1

        # Recur for all the vertices adjacent to this vertex
        for neighbor_id in self.nodes[n_id].neighbors():
            v = node_ids[neighbor_id]
            # If v is not visited yet, then make it a child of u
            # in DFS tree and recur for it
            if disc[v] == -1:
                parent[v] = n_id
                children += 1
                stack.append((n_id, neighbor_id))  # store the edge in stack
                self.bi_cc_rec(
                    neighbor_id,
                    parent,
                    low,
                    disc,
                    stack,
                    node_ids,
                    all_bi_cc,
                    artic_points,
                    disc_time,
                )

                # Check if the subtree rooted with v has a connection to
                # one of the ancestors of u
                # Case 1 -- per Strongly Connected Components Article
                low[u] = min(low[u], low[v])

                # If u is an articulation point, pop
                # all edges from stack until (u, v)
                if parent[u] == -1 and children > 1 or parent[u] != -1 and low[v] >= disc[u]:
                    artic_points.append(n_id)
                    w = -1
                    one_bi_cc = set()
                    while w != (n_id, neighbor_id):
                        w = stack.pop()
                        for n in w:
                            one_bi_cc.add(n)
                    all_bi_cc.append(one_bi_cc)

            elif neighbor_id != parent[u] and low[u] > disc[v]:
                low[u] = min(low[u], disc[v])
                stack.append((n_id, neighbor_id))

    def bicc(self):
        """
        find biconnected components and returns a list of these components in terms of node ids
        """
        disc_time = [0]
        all_bi_cc = list()
        node_ids = dict()
        artic_points = []
        # node_ids_list = [0] * len(self)
        for idx, n_id in enumerate(list(self.nodes.keys())):
            node_ids[n_id] = idx
            # node_ids_list[idx] = n_id

        # Initialize disc and low, and parent arrays
        n_vertices = len(self)
        disc = [-1] * n_vertices
        low = [-1] * n_vertices
        parent = [-1] * n_vertices
        stack = []
        # Call the recursive helper function to
        # find articulation points
        # in DFS tree rooted with vertex 'i'
        for n_id in self.nodes.keys():
            i = node_ids[n_id]
            if disc[i] == -1:
                self.bi_cc_rec(
                    n_id, parent, low, disc, stack, node_ids, all_bi_cc, artic_points, disc_time
                )
            # If stack is not empty, pop all edges from stack
            if stack:
                one_bi_cc = set()
                while stack:
                    w = stack.pop()
                    for n in w:
                        one_bi_cc.add(n)
                all_bi_cc.append(one_bi_cc)
        return all_bi_cc, artic_points
