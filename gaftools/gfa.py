import sys
import logging
from collections import defaultdict

import gzip
import re
import os
from collections import deque
from gaftools.utils import rev_comp, is_correct_tag

E_DIR = {("+", "+"): (1, 0), ("+", "-"): (1, 1), ("-", "+"): (0, 0), ("-", "-"): (0, 1)}


class Node:
    __slots__ = ("id", "seq", "seq_len", "start", "end", "visited", "tags")

    def __init__(self, identifier):
        self.id = identifier
        self.seq = ""
        self.seq_len = 0
        self.start = set()
        self.end = set()
        self.visited = False
        self.tags = dict()

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

    def is_equal_to(self, other, only_topo=False):
        """
        checks if two nodes are equal to each other in terms of edges and other attributes
        """
        all_ats = ["id", "seq", "seq_len", "start", "end", "tags"]
        only_topo_atts = ["id", "start", "end"]
        if only_topo:
            for a in only_topo_atts:
                if not getattr(self, a) == getattr(other, a):
                    return False
        else:
            for a in all_ats:
                if not getattr(self, a) == getattr(other, a):
                    return False
        return True

    def neighbors(self):
        """
        Returns all adjacent nodes' ids to self
        """
        neighbors = [x[0] for x in self.start] + [x[0] for x in self.end]
        return sorted(neighbors)

    def in_direction(self, other, direction):
        """
        returns true if other is a neighbor in that direction, false otherwise
        """
        if direction == 0:
            if other in [x[0] for x in self.start]:
                return True
            return False
        elif direction == 1:
            if other in [x[0] for x in self.end]:
                return True
            return False
        else:
            raise ValueError(
                f"Trying to access a wrong direction in node {self.id}, give 0 for start or 1 for end"
            )

    def children(self, direction):
        """
        returns the children of a node in given direction
        """
        if direction == 0:
            return [x[0] for x in self.start]
        elif direction == 1:
            return [x[0] for x in self.end]
        else:
            raise ValueError(
                f"Trying to access a wrong direction in node {self.id}, give 0 for start or 1 for end"
            )

    def remove_from_start(self, neighbor, side, overlap):
        """
        remove the neighbor edge from the start going to side in neighbor
        """
        assert side in {1, 0}
        try:
            self.start.remove((neighbor, side, overlap))
        except KeyError:
            logging.warning(
                f"Could not remove edge {(neighbor, side, overlap)} from {self.id}'s start as it does not exist"
            )

    def remove_from_end(self, neighbor, side, overlap):
        """
        remove the neighbor edge from the end going to side in neighbor
        """
        assert side in {1, 0}
        try:
            self.end.remove((neighbor, side, overlap))
        except KeyError:
            logging.warning(
                f"Could not remove edge {(neighbor, side, overlap)} from {self.id}'s end as it does not exist"
            )

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

    def to_gfa_line(self, with_seq=True):
        if with_seq:  # in case with_seq was true but there was no seq
            if self.seq == "":
                seq = "*"
            else:
                seq = self.seq
        else:
            seq = "*"
        tags = []
        for tag_id, tag in self.tags.items():
            tags.append(f"{tag_id}:{tag[0]}:{tag[1]}")
        return "\t".join(["S", self.id, seq] + tags)


class GFA:
    """
    Graph object containing the important information about the graph
    """

    __slots__ = ["nodes", "low_memory", "edge_tags", "contigs"]

    def __init__(self, graph_file=None, low_memory=False):
        self.nodes = dict()
        self.edge_tags = dict()
        self.contigs = defaultdict(
            lambda: None
        )  # storing contigs with their SR tag to determine primary contig
        self.low_memory = low_memory
        if graph_file:
            if not os.path.exists(graph_file):
                raise FileNotFoundError
            self.read_graph(gfa_file_path=graph_file, low_memory=low_memory)

    def __len__(self):
        """
        overloading the length function
        """
        return len(self.nodes)

    def __str__(self):
        """
        overloading the string function for printing
        """
        edge_count = 0
        for n in self.nodes.values():
            edge_count += len(n.start)
            edge_count += len(n.end)

        return f"The graph has {len(self)} Nodes and {int(edge_count / 2)} Edges"

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

    def is_equal_to(self, other, only_topo=False):
        """
        a graph equality functions which allows us to compare two graphs together, whether they contain the same information or not

        If only_topo is True, then only the graphs toplogies will be compared, and won't look at sequences and tags of nodes
        """
        if len(self) != len(other):
            return False
        for n_id, node1 in self.nodes.items():
            node2 = other[n_id]
            if node2 is None:  # n_id not in other graph, so not equal graphs
                return False
            if not node1.is_equal_to(node2, only_topo):
                return False
        return True

    def remove_node(self, n_id):
        """
        remove a node and its corresponding edges
        """

        # edges_to_remove = []
        starts = [x for x in self.nodes[n_id].start]
        for n_start in starts:
            overlap = n_start[2]
            # edges_to_remove.append((n_id, 0, n_start[0], n_start[1], overlap))
            self.remove_edge((n_id, 0, n_start[0], n_start[1], overlap))

        ends = [x for x in self.nodes[n_id].end]
        for n_end in ends:
            overlap = n_end[2]
            # edges_to_remove.append((n_id, 1, n_end[0], n_end[1], overlap))
            self.remove_edge((n_id, 1, n_end[0], n_end[1], overlap))

        del self.nodes[n_id]

    def remove_edge(self, edge):
        n1, side1, n2, side2, overlap = edge

        if edge in self.edge_tags:
            del self.edge_tags[edge]

        if side1 == 0:
            self.nodes[n1].remove_from_start(n2, side2, overlap)
        else:
            self.nodes[n1].remove_from_end(n2, side2, overlap)

        if side2 == 0:
            self.nodes[n2].remove_from_start(n1, side1, overlap)
        else:
            self.nodes[n2].remove_from_end(n1, side1, overlap)

    def add_node(self, node_id, seq="", tags=None):
        """
        adds a node to the graph, you need to give at least a node_id
        """
        if not tags:
            tags = []
        node_id = str(node_id)
        if node_id not in self:
            node = Node(node_id)
            node.seq = seq
            node.seq_len = len(seq)
            self[node_id] = node
            # adding the extra tags if any to the node object
            for tag in tags:
                if not is_correct_tag(tag):
                    raise ValueError(
                        f"The tag {tag} for node {node_id} did not match the specifications, check sam specification on tags"
                    )
                tag = tag.split(":")
                # I am adding the tags as key:value, key is tag_name:type and value is the value at the end
                # e.g. SN:i:10 will be {"SN": ('i', '10')}
                self[node_id].tags[tag[0]] = (tag[1], tag[2])  # (type, value)
            if "SN" in self[node_id].tags and "SR" in self[node_id].tags:
                contig_name = self[node_id].tags["SN"][1]
                contig_rank = int(self[node_id].tags["SR"][1])
                if self.contigs[contig_name] is None:
                    self.contigs[contig_name] = contig_rank
                else:
                    assert self.contigs[contig_name] == contig_rank
        else:
            logging.warning(
                f"You are trying to add node {node_id} and it already exists in the graph"
            )

    def add_edge(self, node1, node1_dir, node2, node2_dir, overlap, tags=None):
        assert node1_dir in {"+", "-"}
        assert node2_dir in {"+", "-"}

        node1_dir, node2_dir = E_DIR[(node1_dir, node2_dir)]
        if tags:
            self.edge_tags[(node1, node1_dir, node2, node2_dir)] = tags
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

    def write_graph(
        self, set_of_nodes=None, output_file="output_graph.gfa", append=False, order_bo=False
    ):
        """
        writes a graph file as GFA
        Can be given a set of nodes to only write those nodes with their edges
        if append is true, then it appends to an existing gfa file, otherwise, it overwrites it
        """
        if not output_file.endswith(".gfa"):
            output_file += ".gfa"
        # print("I am here")
        self.write_gfa(
            self,
            set_of_nodes=set_of_nodes,
            output_file=output_file,
            append=append,
            order_bo=order_bo,
        )

    def output_components(self, output_dir):
        """
        writes each connected component in a separate GFA file
        """
        connected_comps = self.all_components()
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
            opened_file = gzip.open(gfa_file_path, "rt")
        elif gfa_file_path.endswith(".gfa"):
            opened_file = open(gfa_file_path, "r")
        else:
            raise ValueError(f"File {gfa_file_path} needs to end with .gfa or .gz")

        for line in opened_file:
            if line.startswith("S"):
                line = line.strip().split("\t")
                assert len(line) >= 3  # must be at least 3 columns for "S id seq"
                if low_memory:
                    self.add_node(line[1], "", line[3:])
                else:
                    self.add_node(line[1], line[2], line[3:])

            elif line.startswith("L"):
                edges.append(line)
        opened_file.close()

        for e in edges:
            e = e.strip().split("\t")
            assert len(e) >= 6  # must be at least 6 columns (L id1 dir1 id2 dir2 overlap)
            e_tags = e[6:]
            e = e[1:6]
            try:
                e[4] = int(e[4][:-1])  # getting overlap
            except ValueError:
                raise ValueError(f"The overlap in {e} has a problem, must be (int)M, e.g. 10M")

            # skip an edge if the node is not there
            if e[0] not in self or e[2] not in self:
                continue
            if not e_tags:  # no tags
                e_tags = [0]
            # TODO this function might be taking too much time when called too many times
            #     reading a graph in GFASubgraph for example is faster
            self.add_edge(*e, e_tags)

    def sort_bo_no(self, set_of_nodes):
        """
        sorts the set_of_nodes based on the BO and NO tags and returns a sorted list
        """
        separate_bubbles = dict()
        for n in set_of_nodes:
            if self[n].tags["BO"][1] not in separate_bubbles:
                separate_bubbles[self[n].tags["BO"][1]] = [n]
            else:
                separate_bubbles[self[n].tags["BO"][1]].append(n)
        bo_ids = []
        for bo, n_list in separate_bubbles.items():
            bo_ids.append(bo)
            separate_bubbles[bo] = sorted(
                separate_bubbles[bo], key=lambda x: int(self.nodes[x].tags["NO"][1])
            )  # sorting the BO bucket by NO

        sorted_set_of_nodes = []
        for bo in sorted(bo_ids):
            for n_id in separate_bubbles[bo]:
                sorted_set_of_nodes.append(n_id)
        return sorted_set_of_nodes

    def write_gfa(
        self, set_of_nodes=None, output_file="output_file.gfa", append=False, order_bo=False
    ):
        """
        Write a gfa out

        :param nodes: Dictionary of nodes object.
        :param set_of_nodes: A list of node ids of the path or nodes we want to generate a GFA file for.
        :param output_file: path to output file
        :param append: if I want to append to a file instead of rewriting it
        """

        if set_of_nodes is None:  # no subgroup of nodes given, then all nodes are outputted
            set_of_nodes = self.nodes.keys()

        if order_bo:
            # sorts first by BO then each BO gets sorted by NO
            sorted_set_of_nodes = self.sort_bo_no(set_of_nodes)
        else:
            sorted_set_of_nodes = set_of_nodes

        if append is False:
            f = open(output_file, "w+")
        else:
            if os.path.exists(output_file):
                f = open(output_file, "a")
            else:
                logging.warning(
                    "Trying to append to a non-existent file\n" "creating an output file"
                )
                f = open(output_file, "w+")

        # going through nodes twice to have the final output as all S line then all L lines
        for n in sorted_set_of_nodes:
            if n not in self.nodes:
                logging.warning("Node {} does not exist in the graph, skipped in output".format(n))
                continue
            line = self.nodes[n].to_gfa_line()
            f.write(line + "\n")

        for n1 in sorted_set_of_nodes:
            if n1 not in self.nodes:
                logging.warning("Node {} does not exist in the graph, skipped in output".format(n1))
                continue

            # writing edges
            edges = []
            for n in self.nodes[n1].start:
                overlap = str(n[2]) + "M"
                if n[0] in set_of_nodes:
                    try:
                        tags = self.edge_tags[(n1, 0, n[0], n[1])]
                    except KeyError:
                        tags = []
                    if tags:
                        if tags[0] == 0:
                            tags = []
                        if n[1] == 0:
                            edge = str(
                                "\t".join(["L", str(n1), "-", str(n[0]), "+", overlap] + tags)
                            )
                            edges.append(edge)
                        else:
                            edge = str(
                                "\t".join(["L", str(n1), "-", str(n[0]), "-", overlap] + tags)
                            )
                            edges.append(edge)

            for n in self.nodes[n1].end:
                overlap = str(n[2]) + "M"

                if n[0] in set_of_nodes:
                    try:
                        tags = self.edge_tags[(n1, 1, n[0], n[1])]
                    except KeyError:
                        tags = []
                    if tags:
                        if tags[0] == 0:
                            tags = []
                        if n[1] == 0:
                            edge = str(
                                "\t".join(["L", str(n1), "+", str(n[0]), "+", overlap] + tags)
                            )
                            edges.append(edge)
                        else:
                            edge = str(
                                "\t".join(["L", str(n1), "+", str(n[0]), "-", overlap] + tags)
                            )
                            edges.append(edge)

            for e in edges:
                f.write(e + "\n")

        f.close()

    def bfs(self, start_id, size=0, reset_visited=True):
        """
        runs BFS from the start node and limits the return to the size given
        if size not given then it keeps going until it runs out of nodes
        """
        if reset_visited:
            self.set_visited(reset_visited)

        if len(self.nodes[start_id].neighbors()) == 0:
            return {start_id}

        if size == 0:
            size = len(self) + 1

        queue = deque()  # deque is faster when popping left
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
        """
        connected_comp = []
        # visited = set()
        for n in self.nodes:
            if not self.nodes[n].visited:
                connected_comp.append(self.find_component(n))
                # visited = visited.union(connected_comp[-1])
        # resetting the visited to false, in case I'm going to use some other algorithm later
        self.set_visited(False)
        return connected_comp

    def set_visited(self, visited=False):
        """
        Resets all nodes to given boolean
        """
        for n in self.nodes.values():
            n.visited = visited

    def graph_from_comp(self, component_nodes):
        new_graph = GFA()
        for n in component_nodes:
            new_node = Node(n)
            new_node.seq = self[n].seq
            new_node.seq_len = self[n].seq_len
            new_node.start = self[n].start
            new_node.end = self[n].end
            new_node.tags = self[n].tags
            new_graph.nodes[n] = new_node
        return new_graph

    def path_exists(self, path):
        """
        Just a sanity check that a path given exists in the graph
        I am assuming that the list of node given as ordered_path taken from the GAF alignment is ordered
        i.e. node 1 parent of node 2, node 2 parent of node 3 and so on
        """
        ordered_path = path.replace(">", ",").replace("<", ",").split(",")
        if ordered_path[0] == "":
            ordered_path = ordered_path[1:]

        for i in range(1, len(ordered_path)):
            current_node = ordered_path[i]
            previous_node = ordered_path[i - 1]
            if current_node in self.nodes[previous_node].neighbors():
                continue
            else:  # some node is not connected to another node in the path
                return False
        return True

    def list_is_path(self, node_list):
        """
        check if a list of nodes create a linear path
        """
        for i in range(1, len(node_list)):
            current_node = node_list[i]
            previous_node = node_list[i - 1]
            if current_node in self.nodes[previous_node].neighbors():
                continue
            else:  # some node is not connected to another node in the path
                return False
        return True

    def return_gfa_path(self, list_of_nodes):
        """
        Given a list of nodes that create a path, returns a string formatted similar to GFA paths
        """

        # I am assuming the list_of_nodes came from get_path, so it's already been checked that it's a correct sorted path
        path = []

        for i in range(len(list_of_nodes) - 1):
            if self.nodes[list_of_nodes[i]].in_direction(list_of_nodes[i + 1], 1):
                path.append(list_of_nodes[i] + "+")
            elif self.nodes[list_of_nodes[i]].in_direction(list_of_nodes[i + 1], 0):
                path.append(list_of_nodes[i] + "-")
            else:
                raise ValueError("Something went wrong in returning a gaf path format")

        if self.nodes[list_of_nodes[-1]].in_direction(list_of_nodes[-2], 0):
            path.append(list_of_nodes[-1] + "+")
        elif self.nodes[list_of_nodes[-1]].in_direction(list_of_nodes[-2], 1):
            path.append(list_of_nodes[-1] + "-")
        else:
            raise ValueError("Something went wrong in returning a gaf path format")

        return ",".join(path)

    def get_path(self, chrom):
        """
        takes a chromosome name (matching the SN tag) and returns the path of that chromosome
        """
        # get all nodes with chrom
        nodes_of_chrom = []
        for n in self.nodes:
            if "SN" in self[n].tags:
                if self[n].tags["SN"][1] == chrom:
                    nodes_of_chrom.append(n)
            else:
                logging.warning(
                    f"Not able to get path for {chrom}, because node {n} doesn't have an SN tag. Stopping! Returning empty list"
                )
                return list()
        # sorting based on SO tags, so we get start to end of chromosome
        sorted_nodes = sorted(nodes_of_chrom, key=lambda x: int(self.nodes[x].tags["SO"][1]))
        # this sorted list should be a path already spelling the chromosome
        if self.list_is_path(sorted_nodes):
            return sorted_nodes
        else:
            logging.warning(
                f"The sorted nodes with SN tag {chrom} did not create a linear path. Stopping! Returning empty list"
            )
            return list()

    def get_contig_length(self, chrom):
        """
        returns the length of the chromosome or contig name
        """
        sorted_nodes = self.get_path(chrom)
        if not sorted_nodes:
            logging.error(
                "Was not able to return the length of the chromosome, check warning message(s)"
            )
            sys.exit(1)
        else:
            return sum([int(self.nodes[x].tags["LN"][1]) for x in sorted_nodes])

    def extract_path(self, path):
        """
        returns the sequences representing that path
        """
        seq = []

        if not self.path_exists(path):
            return ""

        for n in re.findall("[><][^><]+", path):
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

    # def get_bubbles(self, only_simple=False):
    #     """
    #     calls the find_bubbles function and returns the dictionary of all bubbles found
    #     """
    #     from gaftools.find_bubbles import find_bubbles
    #     logging.info("starting to find bubbles in the graph")
    #     bubbles = find_bubbles(self, only_simple)
    #     logging.info(f"finished finding bubbles and found {len(bubbles)}")
    #     return bubbles

    def dfs(self, start_node):
        """
        Performs depth first search from start node given by user
        return the path as a list
        """
        if start_node not in self:
            return []
        if len(self) == 1:
            return [list(self.nodes.keys())[0]]
        if len(self[start_node].neighbors()) == 0:
            return [start_node]

        dfs_out = set()
        ordered_dfs_out = list()
        stack = [start_node]
        while stack:
            s = stack.pop()
            # Note: this is doing membership check on a list, which is O(n) time
            # this is slow, but for now it's ok, might need to change later to a set
            # However, a set is not ordered, so would need to add an ordering function
            if s not in dfs_out:
                dfs_out.add(s)
                ordered_dfs_out.append(s)
            else:
                continue
            for neighbour in self[s].neighbors():
                stack.append(neighbour)
        return ordered_dfs_out

    def biccs(self, set_of_nodes=None):
        """
        This function modelled after Networkx implementation
        https://networkx.org/documentation/networkx-1.9/reference/generated/networkx.algorithms.components.biconnected.biconnected_components.html
        """
        if set_of_nodes is None:
            set_of_nodes = set(self.nodes.keys())

        def edge_stack_to_set(edge_stack):
            """
            turns the edge stack into a set of nodes for the component
            """
            out_set = set()
            for es in edge_stack:
                for n in es:
                    out_set.add(n)
            return out_set

        def next_child(stack_item):
            """
            increments the stack pointer to the next neighbor of stack_item[1] node
            """
            if not stack_item[3]:
                return None
            if stack_item[2] >= len(stack_item[3]):
                return None
            else:
                stack_item[2] += 1
                return stack_item[3][stack_item[2] - 1]

        visited = set()
        components = list()
        artic_points = set()
        for n in set_of_nodes:
            if n in visited:
                continue

            discovery = {n: 0}
            low = {n: 0}
            root_children = 0
            artic_points = set()
            components = []
            visited.add(n)
            edge_stack = []
            edge_stack_loc = dict()
            # stack = [(start, start, iter(G[start]))]
            neighbors = self[n].neighbors()
            stack = [[n, n, 0, neighbors]]
            while stack:
                parent = stack[-1][0]
                child = stack[-1][1]
                nn = next_child(stack[-1])

                if nn:
                    if nn == parent:
                        continue
                    if nn in visited:
                        if discovery[nn] <= discovery[child]:
                            edge_stack.append((child, nn))
                            edge_stack_loc[(child, nn)] = len(edge_stack) - 1
                            low[child] = min(low[child], discovery[nn])
                    else:
                        low[nn] = len(discovery)
                        discovery[nn] = low[nn]
                        visited.add(nn)
                        stack.append([child, nn, 0, self[nn].neighbors()])
                        edge_stack.append((child, nn))
                        edge_stack_loc[(child, nn)] = len(edge_stack) - 1

                elif nn is None:
                    stack.pop()
                    if len(stack) > 1:
                        if low[child] >= discovery[parent]:
                            artic_points.add(parent)
                            cut_point = edge_stack_loc[(parent, child)]
                            # cut_point = edge_stack.index((parent, child))
                            comp = edge_stack_to_set(edge_stack[cut_point:])
                            components.append(comp)
                            del edge_stack[cut_point:]
                            # edge_stack = edge_stack[:cut_point]

                        low[parent] = min(low[parent], low[child])

                    elif stack:
                        root_children += 1
                        cut_point = edge_stack_loc[(parent, child)]
                        comp = edge_stack_to_set(edge_stack[cut_point:])
                        components.append(comp)
                        del edge_stack[cut_point:]

            if root_children > 1:
                artic_points.add(n)

        return components, artic_points
