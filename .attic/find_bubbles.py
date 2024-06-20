"""
This part is taken from BubbleGun, by Fawaz Dabbaghie
Fawaz Dabbaghie, Jana Ebler, Tobias Marschall, BubbleGun: enumerating bubbles and superbubbles in genome graphs, Bioinformatics
https://doi.org/10.1093/bioinformatics/btac448
"""


def find_sb_alg(graph, s, direction, only_simple=False):
    """
    takes the graph and a start node s and add a bubble to the chain
    if one is found if s was the source
    """
    # I tuples of node ids and the direction
    seen = set()
    visited = set()
    nodes_inside = []
    seen.add((s.id, direction))
    # seen.add(s.id)
    S = {(s, direction)}
    while len(S) > 0:
        v = S.pop()
        v = (v[0], v[1])
        visited.add(v[0].id)

        nodes_inside.append(v[0])

        # it's visited so it's not in the seen list anymore
        seen.remove((v[0].id, v[1]))

        # from which direction we are going we take those children
        if v[1] == 0:
            children = v[0].start
        else:
            children = v[0].end

        if len(children) == 0:
            # it's a tip
            break

        for u in children:
            # check where we entered to get children from the other side
            if u[1] == 0:
                u_child_direction = 1
                u_parents = [x[0] for x in graph.nodes[u[0]].start]

            else:
                u_child_direction = 0
                u_parents = [x[0] for x in graph.nodes[u[0]].end]

            if u[0] == s.id:
                # we are in a loop
                S = set()  # so I exit the outer loop too
                break

            # adding child to seen
            # seen.add(u[0])
            if u[1] == 0:
                seen.add((u[0], 1))
            else:
                seen.add((u[0], 0))
            # if all u_parents are visited then we push it into S
            if all(graph.nodes[i].id in visited for i in u_parents):
                S.add((graph.nodes[u[0]], u_child_direction))

        # checking if we finished
        if (len(S) == 1) and (len(seen) == 1):
            t = S.pop()
            nodes_inside.append(t[0])

            if len(nodes_inside) == 2:
                # it's an empty bubble
                # this shouldn't happen if the graph is compacted
                break

            # t[0].visited = True

            # because I'm looking in both directions I end up finding each
            # bubble twice, so I can hash the id of source and sink
            # and see if I've already found it or not
            nodes_inside.remove(s)
            nodes_inside.remove(t[0])
            if only_simple:
                if len(nodes_inside) <= 2:
                    return s, t[0], nodes_inside
                else:
                    return None, None, None
            else:
                return s, t[0], nodes_inside

    return None, None, None


def find_bubbles(graph, only_simple=False):
    """
    main function for finding bubbles
    Takes a graph and fills in the bubble chains
    """

    bubbles = dict()
    for n in graph.nodes.values():
        if n.visited:
            continue
        for d in [0, 1]:  # looking in both direction for each node
            source, sink, inside = find_sb_alg(graph, n, d, only_simple=only_simple)
            if source is not None:
                # if I also skip the source and sink as visited, then I miss another bubble that shares a source or a sink
                # But I can skip the inside of a bubble
                # setting bubble nodes as visited
                # source.visited = True
                # sink.visited = True
                for n in inside:
                    n.visited = True

                if source.id > sink.id:
                    if (source.id, sink.id) not in bubbles:
                        bubbles[(source.id, sink.id)] = [x.id for x in inside]
                else:
                    if (sink.id, source.id) not in bubbles:
                        bubbles[(sink.id, source.id)] = [x.id for x in inside]

    return bubbles
