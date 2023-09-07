"""
Adding a GFA class which will include GFA parsers, graph building, graph functions, graph editing, graph writing and so on

Hopefully used then by all scripts in gaftools

TODO:
1- Now that I have implement add node and add edge, I should change the read graph to start using those functions
2- implement biconnected components so I can get rid of the usage of networkx and keep everything in the GFA class
3- 
"""
import sys
import logging
import gzip
import re
import os
from gaftools.utils import rev_comp, is_correct_tag


class Node:
    __slots__ = ('id', 'seq', 'seq_len', 'start', 'end', 'visited', 'optional')
    def __init__(self, identifier):
        self.id = identifier
        self.seq = ""
        self.seq_len = 0
        self.start = set()
        self.end = set()
        self.visited = False
        self.optional = dict()

    def __sizeof__(self):
        """
        Tries to calculate the size of the node object in bytes
        """
        size = self.id.__sizeof__() + self.seq_len.__sizeof__() + self.visited.__sizeof__()

        if len(self.start) == 0:
            size += self.start.__sizeof__()
        else:
            for i in self.start:
                size += sys.getsizeof(i)

        if len(self.end) == 0:
            size += self.end.__sizeof__()
        else:
            for i in self.end:
                size += sys.getsizeof(i)
        return size

    def __len__(self):
        return self.seq_len

    def neighbors(self):
        """
        Returns all adjacent nodes' ids to self
        """
        neighbors = [x[0] for x in self.start] + [x[0] for x in self.end]
        return sorted(neighbors)

    def in_direction(self, node, direction):
        """
        returns true if node is a neighbor in that direction, false otherwise
        """
        if direction == 0:
            if node in [x[0] for x in self.start]:
                return True
            return False
        elif direction == 1:
            if node in [x[0] for x in self.end]:
                return True
            return False
        else:
            raise ValueError(f"Trying to access a wrong direction in node {self.id}, give 0 for start or 1 for end")

    def children(self, direction):
        """
        returns the children of a node in given direction
        """
        if direction == 0:
            return [x[0] for x in self.start]
        elif direction == 1:
            return [x[0] for x in self.end]
        else:
            raise ValueError(f"Trying to access a wrong direction in node {self.id}, give 0 for start or 1 for end")

    def remove_from_start(self, neighbor, side, overlap):
        """
        remove the neighbor edge from the start going to side in neighbor
        """
        assert side in {1, 0}
        try:
            self.start.remove((neighbor, side, overlap))
        except KeyError:
            logging.warning(f"Could not remove edge {(neighbor, side, overlap)} from {self.id}'s start as it does not exist")

    def remove_from_end(self, neighbor, side, overlap):
        """
        remove the neighbor edge from the end going to side in neighbor
        """
        assert side in {1, 0}
        try:
            self.end.remove((neighbor, side, overlap))
        except KeyError:
            logging.warning(f"Could not remove edge {(neighbor, side, overlap)} from {self.id}'s end as it does not exist")

    def add_from_start(self, neighbor, side, overlap):
        """
        add edge between self.start and neighbor
        """
        assert side in {1, 0}
        self.start.add((neighbor, side, overlap))

    def add_from_end(self, neighbor, side, overlap):
        """
        add edge between self.end and neighbor
        """
        assert side in {1, 0}
        self.end.add((neighbor, side, overlap))

class GFA:
    """
    Graph object containing the important information about the graph
    """
    __slots__ = ['nodes', 'low_memory']

    def __init__(self, graph_file=None, no_seq=False):
        self.nodes = dict()
        self.low_memory = no_seq
        if graph_file:
            if not os.path.exists(graph_file):
                raise FileNotFoundError
        self.read_graph(gfa_file_path=graph_file, low_memory=no_seq)

    def __len__(self):
        """
        overloading the length function
        """
        return len(self.nodes)

    def __str__(self):
        """
        overloading the string function for printing
        """
        return "The graph has {} Nodes".format(len(self.nodes))

    def __contains__(self, key):
        """
        overloading the in operator to check if node exists in graph
        """
        return key in self.nodes

    def __getitem__(self, key):
        """
        overloading the bracket operator
        """
        try:
            return self.nodes[key]
        except KeyError:
            return None

    def __setitem__(self, key, value):
        """
        overloading setting an item in nodes
        """
        if isinstance(value, Node):
            self.nodes[key] = value
        else:
            raise ValueError("the object given to set should be a Node object")

    def __delitem__(self, key):
        """
        overloading deleting item, which removes node and all its edges safely
        """
        self.remove_node(key)

    def reset_visited(self):
        """
        resets all nodes.visited to false
        """
        for n in self.nodes.values():
            n.visited = False

    def remove_node(self, n_id):
        """
        remove a node and its corresponding edges

        TODO: I should stage the edge removals but first make sure it will not fail
        otherwise I'll end up with a node half-removed and will cause a mess
        """
        starts = [x for x in self.nodes[n_id].start]
        for n_start in starts:
            overlap = n_start[2]
            self.remove_edge((n_id, 0, n_start[0], n_start[1], overlap))

        ends = [x for x in self.nodes[n_id].end]
        for n_end in ends:
            overlap = n_end[2]
            self.remove_edge((n_id, 1, n_end[0], n_end[1], overlap))

        del self.nodes[n_id]
        
    def remove_edge(self, edge):
        n1, side1, n2, side2, overlap = edge
        if side1 == 0:
            self.nodes[n1].remove_from_start(n2, side2, overlap)
        else:
            self.nodes[n1].remove_from_end(n2, side2, overlap)

        if side2 == 0:
            self.nodes[n2].remove_from_start(n1, side1, overlap)
        else:
            self.nodes[n2].remove_from_end(n1, side1, overlap)

    def add_node(self, node_id, seq=""):
        if node_id not in self:
            node = Node(node_id)
            node.seq = seq
            node.seq_len = len(seq)
            self[node_id] = node
        else:
            logging.warning(f"You are trying to add node {node_id} and it already exists in the graph")


    def add_edge(self, node1, node1_dir, node2, node2_dir, overlap=0):
        assert node1 in self
        assert node2 in self
        assert node1_dir in {1, 0}
        assert node2_dir in {1, 0}
        if node1_dir == 0:
            self[node1].add_from_start(node2, node2_dir, overlap)
        else:
            self[node1].add_from_end(node2, node2_dir, overlap)

        if node2_dir == 0:
            self[node2].add_from_start(node1, node1_dir, overlap)
        else:
            self[node2].add_from_end(node1, node1_dir, overlap)

    def remove_lonely_nodes(self):
        """
        remove singular nodes with no edges
        """
        nodes_to_remove = [n.id for n in self.nodes.values() if len(n.neighbors()) == 0]
        for i in nodes_to_remove:
            self.remove_node(i)

    def write_graph(self, set_of_nodes=None,
                    output_file="output_graph.gfa",
                    append=False):
        """
        writes a graph file as GFA
        Can be given a set of nodes to only write those nodes with their edges
        if append is true, then it appends to an existing gfa file, otherwise, it overwrites it
        """
        if not output_file.endswith(".gfa"):
            output_file += ".gfa"
        # print("I am here")
        self.write_gfa(self, set_of_nodes=set_of_nodes, output_file=output_file, 
            append=append)

    def output_components(self, output_dir):
        """
        writes each connected component in a separate GFA file
        """
        connected_comps = all_components(self)
        counter = 1
        for cc in connected_comps:
            if len(cc) > 1:
                output_file = output_dir + os.path.sep + "component{}.gfa".format(counter)
                counter += 1
                logging.info("Writing Component {}...".format(output_file))

                self.write_graph(set_of_nodes=cc, output_file=output_file, append=False)


    def read_graph(self, gfa_file_path, low_memory=False):
        """
        Read a gfa file and return a populated graph object
        """
        if not os.path.exists(gfa_file_path):
            logging.error("the gfa file path you gave does not exists, please try again!")
            sys.exit()

        # nodes = dict()
        edges = []
        if gfa_file_path.endswith(".gz"):
            open_file = gzip.open(gfa_file_path, "rt")
        elif gfa_file_path.endswith(".gfa"):
            open_file = open(gfa_file_path, "r")
        else:
            raise ValueError(f"File {gfa_file_path} needs to end with .gfa or .gz")

        for line in open_file:
            if line.startswith("S"):
                line = line.strip().split("\t")
                assert len(line) >= 3  # must be at least "S id seq"
                n_id = str(line[1])
                self.nodes[n_id] = Node(n_id)

                self.nodes[n_id].seq_len = len(line[2])
                if not low_memory:  # don't load the sequence
                    self.nodes[n_id].seq = str(line[2]).strip()
                # adding the extra tags if any to the node object
                if len(line) > 3:
                    for tag in line[3:]:
                        if not is_correct_tag(tag):
                            raise ValueError(f"The tag {tag} did not match the specifications, check sam specification on tags")
                        tag = tag.split(":")
                        # I am adding the tags as key:value, key is tag_name:type and value is the value at the end
                        self[n_id].optional[f"{tag[0]}:{tag[1]}"] = tag[2]
                        # self.nodes[n_id].optional.append(tag)
                    # self.nodes[n_id].optional = "\t".join(line[3:])

            elif line.startswith("L"):
                edges.append(line)
        open_file.close()

        for e in edges:
            line = e.split()
            k = line[1]
            if k not in self.nodes:  # if the edge is there but not the node
                continue
            try:
                overlap = int(line[5][:-1])
            except:
                raise ValueError(f"Overlap should be represented as (int)M, e.g. 10M, could not recognize {overlap}")
            neighbor = line[3]
            if neighbor not in self.nodes:
                continue

            from_start = False
            to_end = False
            if line[2] == "-":
                from_start = True
            if line[4] == "-":
                to_end = True

            if from_start and to_end:  # from start to end L x - y -
                if (neighbor, 1, overlap) not in self.nodes[k].start:
                    self.nodes[k].start.add((neighbor, 1, overlap))
                if (k, 0, overlap) not in self.nodes[neighbor].end:
                    self.nodes[neighbor].end.add((k, 0, overlap))

            elif from_start and not to_end:  # from start to start L x - y +

                if (neighbor, 0, overlap) not in self.nodes[k].start:
                    self.nodes[k].start.add((neighbor, 0, overlap))

                if (k, 0, overlap) not in self.nodes[neighbor].start:
                    self.nodes[neighbor].start.add((k, 0, overlap))

            elif not from_start and not to_end:  # from end to start L x + y +
                if (neighbor, 0, overlap) not in self.nodes[k].end:
                    self.nodes[k].end.add((neighbor, 0, overlap))

                if (k, 1, overlap) not in self.nodes[neighbor].start:
                    self.nodes[neighbor].start.add((k, 1, overlap))

            elif not from_start and to_end:  # from end to end L x + y -
                if (neighbor, 1, overlap) not in self.nodes[k].end:
                    self.nodes[k].end.add((neighbor, 1, overlap))

                if (k, 1, overlap) not in self.nodes[neighbor].end:
                    self.nodes[neighbor].end.add((k, 1, overlap))

    def write_gfa(self, set_of_nodes=None, output_file="output_file.gfa", append=False):
        """
        Write a gfa out

        :param nodes: Dictionary of nodes object.
        :param set_of_nodes: A list of node ids of the path or nodes we want to generate a GFA file for.
        :param output_file: path to output file
        :param append: if I want to append to a file instead of rewriting it
        """
        # nodes = self.nodes

        if set_of_nodes is None:
            set_of_nodes = self.nodes.keys()

        if append is False:
            f = open(output_file, "w+")
        else:
            if os.path.exists(output_file):
                f = open(output_file, "a")
            else:
                logging.warning("Trying to append to a non-existent file\n"
                                "creating an output file")
                f = open(output_file, "w+")
        for n1 in set_of_nodes:
            if n1 not in self.nodes:
                logging.warning("Node {} does not exist in the graph, skipped in output".format(n1))
                continue

            if self.nodes[n1].optional:
                tags = []
                for tag_name, value in self.nodes[n1].optional.items():
                    tags.append(f"{tag_name}:{value}")
                line = str("\t".join(["S", str(n1), self.nodes[n1].seq, "\t".join(tags)]))
                # line = str("\t".join(("S", str(n1), nodes[n1].seq, nodes[n1].optional)))
            else:
                line = str("\t".join(["S", str(n1), self.nodes[n1].seq + "\n"]))

            f.write(line+"\n")

            # writing edges
            edges = []
            # overlap = str(graph.k - 1) + "M\n"

            for n in self.nodes[n1].start:
                overlap = str(n[2]) + "M\n"

                if n[0] in set_of_nodes:
                    if n[1] == 0:
                        edge = str("\t".join(("L", str(n1), "-", str(n[0]), "+", overlap)))
                        edges.append(edge)
                    else:
                        edge = str("\t".join(("L", str(n1), "-", str(n[0]), "-", overlap)))
                        edges.append(edge)

            for n in self.nodes[n1].end:
                overlap = str(n[2]) + "M\n"

                if n[0] in set_of_nodes:
                    if n[1] == 0:
                        edge = str("\t".join(("L", str(n1), "+", str(n[0]), "+", overlap)))
                        edges.append(edge)
                    else:
                        edge = str("\t".join(("L", str(n1), "+", str(n[0]), "-", overlap)))
                        edges.append(edge)

            for e in edges:
                f.write(e)

        f.close()

    def bfs(self, start_id, size, reset_visited=True):
        if reset_visited:
            self.reset_visited()

        if len(self.nodes[start_id].neighbors()) == 0:
            return {start_id}

        queue = deque()
        queue.append(start_id)
        self.nodes[start_id].visited = True
        neighborhood = set()
        while (len(neighborhood) <= size) and len(queue) > 0:
            s = queue.popleft()
            neighborhood.add(s)
            self.nodes[s].visited = True
            for n in self.nodes[s].neighbors():
                if not self.nodes[n].visited:
                    queue.append(n)
                    self.nodes[n].visited = True
        return neighborhood

    def find_component(self, start_node):
        """
        Find component in the graph starting from given node

        :param graph: is a graph object from class Graph
        :param start_node: The start node of BFS search.
        :return: a list of node ids for the component
        """
        queue = []
        cc = set()

        # visited = set()
        queue.append(start_node)
        self.nodes[start_node].visited = True
        # visited.add(start_node)
        neighbors = self.nodes[start_node].neighbors()

        if len(neighbors) == 0:
            cc.add(start_node)
            return cc

        while len(queue) > 0:
            start = queue.pop()
            if start not in cc:
                cc.add(start)
            else:
                continue

            # visited.add(start)
            self.nodes[start].visited = True
            neighbors = self.nodes[start].neighbors()
            for n in neighbors:
                if not self.nodes[n].visited:
                    queue.append(n)

        return cc

    def all_components(self):
        """
        find all connected components in the graph

        :params graph: is a graph object from class Graph
        :return: list of set of components
        """
        connected_comp = []
        # visited = set()
        for n in self.nodes:
            if not graph.nodes[n].visited:
                connected_comp.append(find_component(self, n))
                # visited = visited.union(connected_comp[-1])
        return connected_comp

    def reset_visited(self, visited=False):
        """
        Resets all nodes to given boolean
        """
        for n in self.nodes.values():
            n.visited = visited

    def path_exists(self, path):
        """
        Just a sanity check that a path given exists in the graph
        I am assuming that the list of node given as ordered_path taken from the GAF alignment is ordered
        i.e. node 1 parent of node 2, node 2 parent of node 3 and so on
        """
        print("wkjdf")

        # might be hacky, will think of a better way later
        ordered_path = path.replace(">", ",").replace("<", ",").split(",")
        if ordered_path[0] == "":
            ordered_path = ordered_path[1:]

        for i in range(1, len(ordered_path)):
            current_node = ordered_path[i]
            previous_node = ordered_path[i-1]
            if current_node in self.nodes[previous_node].neighbors():
                continue
            else:
                logging.error(f"in path {ordered_path}, node {current_node} is not a neighbor of {previous_node}, "
                              "no continuous walk through the path")
                # print("The path is not a valid path")
                return False
        return True

    def extract_path(self, path):
        """
        returns the sequences representing that path
        """
        seq = []

        if not self.path_exists(path):
            return ""

        for n in re.findall('[><][^><]+', path):
            if n[1:] not in self:
                logging.error(f"The node {n[1:]} in path {path} does not seem to exist in this GFA")
                return ""
            # print(n)
            if n.startswith(">"):
                seq.append(self.nodes[n[1:]].seq)
            elif n.startswith("<"):
                seq.append(rev_comp(self.nodes[n[1:]].seq))
                # seq.append("".join([reverse_complement[x] for x in self.nodes[n[1:]].seq[::-1]]))
            else:
                logging.error(f"Some error happened where a node {n} doesn't start with > or <")
                return ""

        return "".join(seq)
