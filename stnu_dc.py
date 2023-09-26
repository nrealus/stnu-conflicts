"""
The functions defined in this module allow us to check dynamic
controllability of a STNU and, if it is uncontrollable, to extract
"conflicts" which provide a reason to why it is dynamically uncontrollable.

The algorithm for dynamic controllability checking and conflicts
extraction is based on [1], which in turn is based on Morris' O(n^3)
dynamic controllability checking algorithm [2].

[1]: N.Bhargava et al. *Faster Conflict Generation for Dynamic Controllability* (2017)
[2]: P.Morris *Dynamic Controllability and Dispatchability Relationships* (2014)
"""

from __future__ import annotations

#################################################################################

from heapdict import heapdict
from typing import Dict, List, NamedTuple, Optional, Set, Tuple, Iterable
from enum import Enum

#################################################################################
#
#################################################################################

Node = str
"""A type representing the nodes in an STNU (i.e. timepoints)."""

Weight = float
"""
A type representing the weight of an edge in an STNU
(i.e. "distance" between two timepoints).
"""

class Label(NamedTuple):
    """
    A type representing the label of an edge in an STNU.
    Labels with `kind = Kind.NONE` must have `node = None`.
    Similarly, nodes with `kind = Kind.LOWERCASE` or
    `kind = Kind.UPPERCASE` must have `node != None`
    """
    class Kind(Enum):
        NONE = 0
        UPPERCASE = 1
        LOWERCASE = 2
    kind: Kind
    node: Optional[Node]

    @classmethod
    def nocase(cls):
        return Label(Label.Kind.NONE, None)

    @classmethod
    def lowercase(cls,
        node: Node
    ):
        return Label(Label.Kind.LOWERCASE, node)

    @classmethod
    def uppercase(cls,
        node: Node
    ):
        return Label(Label.Kind.UPPERCASE, node)

InverseLDG = Dict[Node, Dict[Node, Tuple[Weight, Label]]]
"""
Inverse (or reverse) "view" of a "Labeled Distance Graph" (LDG).

The key of the outer dictionary represents the target (not source!)
node of an STNU edge. The inner dictionary represents the source as
well as the weight and label of this STNU edge.

Example: `inverse_ldg[tgt_node] = { src_node: (w1, l1), src2: (w2, l2), ...}`
"""

Edge = Tuple[Node, Node, Weight, Label]
"""
Representation of an STNU edge as a tuple:

source node, target node, weight, label
"""
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class STNU(NamedTuple):
    inverse_ldg: InverseLDG
    """
    Contains the (inverse) labeled distance graph
    representation of the STNU.
    """
    contingent_links: Set[Tuple[Node, Node]]
    """
    Contains the contingent links of the STNU.
    """

    @classmethod
    def from_links(cls,
        ordinary: List[Tuple[Node, Weight | Tuple[Weight, Weight], Node]],
        contingent: List[Tuple[Node, Tuple[Weight, Weight], Node]],
    ):

        stnu = STNU({}, { (src_node, tgt_node) 
                    for (src_node, (_, _), tgt_node) in contingent})

        for (src_node, (l, u), tgt_node) in contingent:
            
            if ((tgt_node in stnu.inverse_ldg and src_node in stnu.inverse_ldg[tgt_node])
                or (src_node in stnu.inverse_ldg and tgt_node in stnu.inverse_ldg[src_node])
            ):
                raise ValueError("Input edges not unique: there is more than one edge defined between two nodes.")

            stnu.inverse_ldg.setdefault(tgt_node, {})[src_node] = (l, Label.lowercase(tgt_node))
            stnu.inverse_ldg.setdefault(src_node, {})[tgt_node] = (-u, Label.uppercase(tgt_node))

        for (src_node, ws, tgt_node) in ordinary:

            if ((tgt_node in stnu.inverse_ldg and src_node in stnu.inverse_ldg[tgt_node])
                or (src_node in stnu.inverse_ldg and tgt_node in stnu.inverse_ldg[src_node])
            ):
                raise ValueError("Input edges not unique: there is more than one edge defined between two nodes.")

            if isinstance(ws, tuple):
                l, u = ws
                stnu.inverse_ldg.setdefault(tgt_node, {})[src_node] = (u, Label.nocase())
                stnu.inverse_ldg.setdefault(src_node, {})[tgt_node] = (-l, Label.nocase())
            else:
                stnu.inverse_ldg.setdefault(tgt_node, {})[src_node] = (ws, Label.nocase())

        return stnu

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def convert_stnu_to_normal_form(
    stnu: STNU,
) -> Tuple[STNU, Dict[Node, Tuple[Node, Node]]]:
    """
    Converts the given STNU into a new STNU in "normal form".

    The normal form of an original STNU, is a STNU where the original
    contingent links have been split in two links: one ordinary (or
    requirement) link, and one contingent link. Both of these links
    are connected using a newly added helper node. The length of the
    newly created ordinary link is fixed / constant, and the lower
    bound of the duration of the contingent link is 0.

    Example: The contingent link A ==[l,u]==> C is transformed into
    these two links: A --[l,l]--> A'' ==[0,u-l]==> C.

    Returns:
        - A new STNU, corresponding to the "normal form" of the original STNU.
        - A mapping from the helper nodes in the new STNU to the
        original STNU contingent links that these nodes helped to split.

    """
    
    normal_form_stnu = STNU({ node : stnu.inverse_ldg[node].copy()
                              for node in stnu.inverse_ldg },
                            set())

    normal_form_stnu_helper_nodes: Dict[Node, Tuple[Node, Node]] = {}

    helpers_counter: int = 0

    for src_node, tgt_node in stnu.contingent_links:

        weight_src_to_tgt = stnu.inverse_ldg[tgt_node][src_node][0]
        weight_tgt_to_src = -stnu.inverse_ldg[src_node][tgt_node][0]

        normal_form_stnu.inverse_ldg[tgt_node].pop(src_node)
        normal_form_stnu.inverse_ldg[src_node].pop(tgt_node)

        helpers_counter += 1
        helper_node = Node(src_node+"''"+str(helpers_counter))

        normal_form_stnu.inverse_ldg[helper_node] = {}

        normal_form_stnu.inverse_ldg[helper_node][src_node] = (weight_src_to_tgt, Label.nocase())
        normal_form_stnu.inverse_ldg[src_node][helper_node] = (-weight_src_to_tgt, Label.nocase())
        normal_form_stnu.inverse_ldg[tgt_node][helper_node] = (Weight(0), Label.lowercase(tgt_node))
        normal_form_stnu.inverse_ldg[helper_node][tgt_node] = (weight_src_to_tgt - weight_tgt_to_src,
                                                               Label.uppercase(tgt_node))

        normal_form_stnu.contingent_links.add((helper_node, tgt_node))
        normal_form_stnu_helper_nodes[helper_node] = (src_node, tgt_node)

    return (normal_form_stnu, normal_form_stnu_helper_nodes)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def check_stnu_dc(
    stnu: STNU,
) -> Tuple[bool, List[List]]:
    """
    Checks whether the given STNU is dynamically controllable (DC),
    and if not, provides a reason to why it isn't, in the form of a
    set of "conflicts" (which are specific paths from the semi-reducible
    negative cycle found to prove that the STNU is not DC).

    Implements the top level of Morris' O(n^3) DC-checking algorithm,
    with slight modifications as in Bhargava et al. (Algorithm 1).

    Returns:
        - Whether the STNU is dynamically controllable
        - If the STNU is not dynamically controllable, a set of "conflicts".
    """
    
    # The inverse labeled distance graph of the normal
    # form STNU may be modified by `dc_disjktra`. Newly added edges
    normal_form_stnu, normal_form_stnu_helper_nodes = convert_stnu_to_normal_form(stnu)

    _normal_form_stnu_inverse_ldg_copy = { node: incoming.copy()
                                           for node, incoming in normal_form_stnu.inverse_ldg.items() }

    # The booleans in `negative_nodes` indicate whether a negative node
    # (i.e. a node with a negative edge directed at it) has already been
    # (successfully) processed during `dc_dijkstra`, meaning it that there
    # are no semi-reducible negative cycle starting (backwards) from that node.
    negative_nodes: Dict[Node, bool] = { node: False
                                          for node, incoming in normal_form_stnu.inverse_ldg.items()
                                          for _, (weight, _) in incoming.items() if weight < 0 }

    # The edges newly added to the (normal form) STNU during runs of `dc_dijkstra`.
    # Note that updated already existing edges won't be recorded.
    novel_edges: Set[Tuple[Node, Node]] = set()

    # `distances_to` stores the shortest distance from a source node
    # (key of the inner dictionary) to a target node (key of the outer 
    # dictionary) as well as the first edge of the shortest path between them.
    distances_to: Dict[Node, Dict[Node, Tuple[Weight, Optional[Edge]]]] = {}

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    for start_node, was_processed in negative_nodes.items():
        
        if was_processed:
            continue

        dc, raw_paths, _ = dc_dijkstra(start_node,
                                       normal_form_stnu.inverse_ldg,
                                       novel_edges,
                                       negative_nodes,
                                       distances_to,
                                       [])

        if not dc:

            # The paths in `raw_paths` form a semi-reducible negative cycle.
            # These paths result from the application of a series of
            # "reduction rules" (see STNU theory) on lowercase edges.
            #
            # The "conflicts" are the full semi-reducible negative cycle
            # itself (or rather the path of edges composing it), as well
            # as these paths / sets of edges described above,
            # without their first edge.

            # But first, before we extract the conflicts, we need to
            # process these raw paths. This is done in two steps, for
            # each raw path:
            #
            # - Decomposing those of their edges that were
            #   not present in the initial STNU in normal form.
            # 
            # - Converting the edges from the initial normal form STNU
            #   to edges of the original, not necessarily normal form STNU.

            processed_paths: List[List[Tuple[Node, Node]]] = []

            for raw_path in raw_paths:
                processed_paths.append(process_raw_path(raw_path,
                                                        novel_edges,
                                                        distances_to,
                                                        normal_form_stnu_helper_nodes))
            
            conflicts: List[List[Tuple[Node, Node]]] = [[]]

            for processed_path in processed_paths:
                conflicts[0].extend(processed_path)

            if len(processed_paths) > 1:
                for processed_path in processed_paths:
                    if len(processed_path) > 1:
                        conflicts.append(processed_path[1:])

            print("----")
            print(raw_paths)

            print("--")
            print(processed_paths)

            print("--")
            print(conflicts)
            
            ws = []
            for conflict in conflicts:
                w = 0
                for src_node, tgt_node in conflict:
                    w += stnu.inverse_ldg[tgt_node][src_node][0]
                ws.append(w)
            print(ws)
            
            print("----")

            return False, conflicts

    return True, []

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def dc_dijkstra(
    start_node: Node,
    normal_form_stnu_inverse_ldg: InverseLDG,
    novel_edges: Set[Tuple[Node, Node]],
    negative_nodes: Dict[Node, bool],
    distances_to: Dict[Node, Dict[Node, Tuple[Weight, Optional[Edge]]]],
    call_stack: List[Node],
) -> Tuple[bool,
           List[List[Tuple[Node, Node]]],
           Optional[Node]]:
    """
    The main function of the DC-checking algorithm.

    Implements Bhargava et al.'s Algorithm 2, which is based on
    Morris' determineDC procedure.
    """

    call_stack.append(start_node)

    if start_node in call_stack[:-1]:
        call_stack.pop()
        return False, [], start_node

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Items type: Tuple[Node, Label]
    # Priority type: Weight
    priority_queue = heapdict()

    distances_to_start_node_from: Dict[Node, Tuple[Weight, Optional[Edge]]] = \
        { start_node : (Weight(0), None) }

    for e_source_node, (e_weight, e_label) in normal_form_stnu_inverse_ldg[start_node].items():
        edge: Edge = (e_source_node, start_node, e_weight, e_label)

        if e_weight < 0:
            distances_to_start_node_from[e_source_node] = (e_weight, edge)
            priority_queue[(e_source_node, e_label)] = e_weight

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    distances_to[start_node] = distances_to_start_node_from

    while priority_queue:

        (node, label), weight = priority_queue.popitem()

        node: Node
        label: Label
        weight: Weight

        if weight >= 0:
            if (start_node not in normal_form_stnu_inverse_ldg
                or node not in normal_form_stnu_inverse_ldg[start_node]
            ):
                novel_edges.add((node, start_node))
            normal_form_stnu_inverse_ldg[start_node][node] = (weight, Label.nocase())
            continue

        if node in negative_nodes:

            srnc_found, raw_paths, end_node = dc_dijkstra(node,
                                                          normal_form_stnu_inverse_ldg,
                                                          novel_edges,
                                                          negative_nodes,
                                                          distances_to,
                                                          call_stack)
            if not srnc_found:

                if end_node is not None:

                    raw_path: List[Tuple[Node, Node]] = build_raw_path(node,
                                                                       start_node,
                                                                       distances_to)
                    raw_paths.append(raw_path)

                if end_node == start_node:
                    end_node = None

                call_stack.pop()
                return False, raw_paths, end_node

        for e_source_node, (e_weight, e_label) in normal_form_stnu_inverse_ldg[node].items():
            edge: Edge = (e_source_node, node, e_weight, e_label)

            if e_weight < 0:
                continue

            if (e_label.kind == Label.Kind.LOWERCASE
                and label.kind != Label.Kind.NONE
                and e_label.node == label.node 
            ):
                continue

            w = weight + e_weight
            
            if (e_source_node not in distances_to_start_node_from
                or w < distances_to_start_node_from[e_source_node][0]
            ):
                distances_to_start_node_from[e_source_node] = (w, edge)
                priority_queue[(e_source_node, label)] = w

    negative_nodes[start_node] = True

    call_stack.pop()
    return True, [], None

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def build_raw_path(
    src_node: Node,
    tgt_node: Node,
    distances_to: Dict[Node, Dict[Node, Tuple[Weight, Optional[Edge]]]],
) -> List[Tuple[Node, Node]]:
    """
    Builds the shortest path from `src_node` to `tgt_node`,
    by "hopping" from node to node using `distances_to`.

    Note that the path may go through edges absent from
    the initial normal form STNU. These edges would
    be added during the run of `dc_dijkstra`.
    """

    path = []
    cur_node = src_node

    while True:
        _, edge = distances_to[tgt_node][cur_node]
        assert edge is not None

        next_node = edge[1]

        path.append((cur_node, next_node))

        if next_node == tgt_node:
            break
        cur_node = next_node

    return path

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def process_raw_path(
    raw_path: List[Tuple[Node, Node]],
    novel_edges: Set[Tuple[Node, Node]],
    distances_to: Dict[Node, Dict[Node, Tuple[Weight, Optional[Edge]]]],
    normal_form_stnu_helper_nodes: Dict[Node, Tuple[Node, Node]],
) -> List[Tuple[Node, Node]]:
    """
    Processes a raw path in two passes:

    1. Decomposing the edges that were not present
    in the initial normal form STNU

    2. Converting the edges (that now all match to edges in
    the normal form STNU) to edges of the original, not
    necessarily normal form STNU. In other words, remove
    helper nodes (however it is slightly more complicated
    than just removing them! see comment below.)
    """

    normal_form_edges_path = []

    for src_node, tgt_node in raw_path:

        stack = [(src_node, tgt_node)]

        while stack:
            src_node, tgt_node = stack.pop()

            if (src_node,tgt_node) not in novel_edges:
                normal_form_edges_path.append((src_node, tgt_node))
            
            else:

                ext = []
                _, e = distances_to[tgt_node][src_node]
                assert e is not None

                while True:
                    ext.append((e[0], e[1]))

                    if e[1] == tgt_node:
                        break

                    _, e = distances_to[tgt_node][e[1]]
                    assert e is not None
                
                stack.extend(reversed(ext))

    # To remove helper nodes in the path, the following rules are applied:
    # 
    # - If the first edge of the path starts from a helper node,
    # replace that helper node with the source node of the original
    # contingent link which was split using that helper node.
    # 
    # - Similarly, if the last edge of the path ends with a helper node,
    # replace it with the target node of the original contingent link
    # which was split using that helper node.
    # 
    # - If any (except the last one) edge of the path ends with a helper node,
    # check whether this edge's source node is the same as the next edge's
    # target node (i.e. whether these two edges are "symmetric"). If they
    # are symmetric, we can "unify" these two edges into one, "skipping"
    # the helper node inbetween. If they are symmetric, we extend them both
    # (symmetrically) up to the target node of the original contingent link.

#    print(normal_form_edges_path)

    original_edges_path = []

    # TODO: The following works but... hic sunt dracones...!
    # should come back to this at some point and try to simplify / clarify this.

    n = len(normal_form_edges_path)
    for i, (src_node, tgt_node) in enumerate(normal_form_edges_path):

        stack = [(src_node, tgt_node)]

        if src_node in normal_form_stnu_helper_nodes:

            if i == 0:
                s, t = normal_form_stnu_helper_nodes[src_node]
                if s != tgt_node:
                    original_edges_path.append((s, tgt_node))
                else:
                    original_edges_path.append((t, tgt_node))

        elif tgt_node in normal_form_stnu_helper_nodes:

            if i == n-1:
                s, _ = normal_form_stnu_helper_nodes[tgt_node]
                original_edges_path.append((src_node, s))

            elif i < n-1:
                
                s, t = normal_form_stnu_helper_nodes[tgt_node]
                ss, tt = normal_form_edges_path[i+1]
                
                if src_node == tt:
                    if src_node == s:
                        original_edges_path.append((src_node, t))
                        original_edges_path.append((t, src_node))
                    else:
                        original_edges_path.append((src_node, s))
                        original_edges_path.append((s, src_node))
                else:
                    original_edges_path.append((src_node, tt))

                continue

            else:
                assert False

        else:
            original_edges_path.append((src_node, tgt_node))

    return original_edges_path

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import unittest

class Tests(unittest.TestCase):

    def test01_non_dc(self):
        stnu = STNU.from_links([("C", -2 ,"A"),
                                ("A", 3, "B"),
                                ("B", -2, "C")],
                               [])
        results = check_stnu_dc(stnu)
        self.assertFalse(results[0])

    def test02_dc(self):
        stnu = STNU.from_links([("C", -1 ,"A"),
                                ("A", 3, "B"),
                                ("B", -2, "C")],
                               [])
        results = check_stnu_dc(stnu)
        self.assertTrue(results[0])

    def test03_dc(self):
        stnu = STNU.from_links([("A", 2, "C"),
                                ("C", 2, "B")],
                               [("A", (0, 3), "B")])
        results = check_stnu_dc(stnu)
        self.assertTrue(results[0])

    def test04_non_dc(self):
        stnu = STNU.from_links([("E", (0.5, 1), "H"),
                                ("H", (0, 2), "I"),
                                ("G", (0, 3), "I")],
                               [("D", (6, 11), "E"),
                                ("F", (2, 4), "G")])
        results = check_stnu_dc(stnu)
        self.assertFalse(results[0])

    def test05_dc(self):
        stnu = STNU.from_links([("X", (7, 12), "C0"),
                                ("C1", (1, 11), "C0")],
                               [("A0", (1, 3), "C0"),
                                ("A1", (1, 10), "C1")])
        results = check_stnu_dc(stnu)
        self.assertTrue(results[0])

    def test06_dc(self):
        stnu = STNU.from_links([("X", (29, 48), "C0"),
                                ("C1", (1, 8), "C0"),
                                ("C2", (7, 37), "C0")],
                               [("A0", (1, 3), "C0"),
                                ("A1", (1, 10), "C1"),
                                ("A2", (1, 36), "C2")])
        results = check_stnu_dc(stnu)
        self.assertTrue(results[0])

    def test07_non_dc(self):
        stnu = STNU.from_links([("E1", (0, 50), "E3"),
                                ("E2", (40, 45), "E3")],
                               [("E1", (20, 30), "E2")])
        results = check_stnu_dc(stnu)
        self.assertFalse(results[0])

    def test08_non_dc(self):
        stnu = STNU.from_links([("E2", (1, 2), "E3"),
                                ("E1", (1, 10), "E2")],
                               [("E1", (1, 8), "E3")])
        results = check_stnu_dc(stnu)
        self.assertFalse(results[0])

    def test09_dc(self):
        stnu = STNU.from_links([("10", (399, 650), "1"),
                                ("9", (109, 468), "1"),
                                ("7", (29, 128), "1"),
                                ("5", (7, 34), "1"),
                                ("3", (1, 8), "1")],
                               [("0", (1, 3), "1"),
                                ("2", (1, 10), "3"),
                                ("4", (1, 36), "5"),
                                ("6", (1, 130), "7"),
                                ("8", (1, 470), "9")])
        results = check_stnu_dc(stnu)
        self.assertTrue(results[0])

    def test10_non_dc(self):
        stnu = STNU.from_links([("E1", (1, 100), "E2"),
                                ("E2", (0, 100), "E5"),
                                ("E2", (1, 100), "E3"),
                                ("E3", (1.5, 100), "E4"),
                                ("E1", (0, 3.5), "E4")],
                               [("E1", (0.6294, 18.8554), "E5")])
        results = check_stnu_dc(stnu)
        self.assertFalse(results[0])

    def test11_non_dc(self):
        stnu = STNU.from_links([("E4", (-1, 3), "E2"),
                                ("E5", (2, 4), "E2")],
                               [("E1", (0, 2), "E2"),
                                ("E3", (0, 3), "E4")])
        results = check_stnu_dc(stnu)
        self.assertFalse(results[0])

    def test12_non_dc(self):
        stnu = STNU.from_links([("E2", (0, 0), "E4"),
                                ("E3", (0, 0), "E4")],
                               [("E1", (0, 1), "E2"),
                                ("E1", (0, 1), "E3")])
        results = check_stnu_dc(stnu)
        self.assertFalse(results[0])

    def test13_dc(self):
        stnu = STNU.from_links([("E1", (0, 2), "E2"),
                                ("E2", (0, 2), "E3")],
                               [("E1", (0, 4), "E3")])
        results = check_stnu_dc(stnu)
        self.assertTrue(results[0])

    def test14_non_dc(self):
        stnu = STNU.from_links([("E3", (0, 2), "E2")],
                               [("E1", (0, 10), "E2"),
                                ("E1", (0, 8), "E3")])
        results = check_stnu_dc(stnu)
        self.assertFalse(results[0])

    def test15_non_dc(self):
        stnu = STNU.from_links([("E2", (-1, 1000000), "E3"),
                                ("E4", (0, 1000000), "E5"),
                                ("E6", (0, 1000000), "E7"),
                                ("E2", (5, 18), "E7")],
                               [("E1", (3, 1000000), "E2"),
                                ("E3", (1, 5.5), "E4"),
                                ("E5", (10, 14.5), "E6")])
        results = check_stnu_dc(stnu)
        self.assertFalse(results[0])

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == "__main__":
    unittest.main()
