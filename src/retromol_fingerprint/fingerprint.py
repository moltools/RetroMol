"""Module containing fingerprinting functionality of module graphs."""

import hashlib
import math
from collections import Counter
from typing import Any

import numpy as np


try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False


def module_graph_fingerprint(
    graph: "nx.Graph",
    substitution_matrix: dict[str, dict[str, float]],
    n_bits: int = 2048,
    radius: int = 2,
    node_name_attr: str = "name",
    # embedding_dim: int = 16,
    embedding_dim: int = 4,
    count_similar_edges: bool = True,
    counted: bool = False,
) -> np.ndarray:
    """
    Compute a fixed-length fingerprint of a module graph, using a custom 
    algorithm that incorporates substitution-aware module embeddings derived 
    from the provided matrix, as well as structural features of the graph.

    :param graph: The module graph to fingerprint. Nodes must have a string attribute
        specified by `node_name_attr` that indicates the module identity.
    :param substitution_matrix: A nested dictionary mapping module identities to their
        pairwise substitution scores. Must include all module identities present in the graph.
    :param n_bits: The length of the output fingerprint in bits.
    :param radius: The number of iterations of neighborhood aggregation to perform,
        similar to a Weisfeiler-Lehman graph kernel. Higher values capture more global structure
        but are more expensive to compute.
    :param node_name_attr: The node attribute key to use for looking up module identities.
    :param embedding_dim: The dimensionality of the substitution-aware module embeddings to derive
        from the substitution matrix. Should be less than n_bits for effective hashing.
    :param count_similar_edges: Whether to hash features that count similar edges (same module
        pairs) together, or just hash presence/absence of edge types. Counting can help distinguish
        between graphs with different numbers of similar edges, but may also make the fingerprint more sensitive to
        small changes.
    :return: A fixed-length numpy array representing the fingerprint of the graph. The dtype is uint16 if `counted` is True
        to allow for counting features, otherwise uint8 for a binary fingerprint.
    :raise ValueError: If the graph contains nodes missing the required name attribute, or if the substitution matrix is 
        missing entries for any module names found in the graph, or if embeddings cannot be derived from the
        substitution matrix.
    :raise ImportError: If NetworkX is not available, since it is required for graph processing.
    """
    if not NETWORKX_AVAILABLE:
        raise ImportError("NetworkX is required for module graph fingerprinting! Please install it with `pip install networkx`!")

    def hash_to_index(*parts: Any) -> int:
        text = "|".join(map(str, parts))
        digest = hashlib.blake2b(text.encode("utf-8"), digest_size=16).digest()
        return int.from_bytes(digest, "big") % n_bits

    def quantize(vec: np.ndarray, decimals: int = 2) -> tuple[float, ...]:
        return tuple(np.round(vec.astype(float), decimals=decimals))
    
    if graph.number_of_nodes() == 0:
        return np.zeros(n_bits, dtype=np.uint16 if counted else np.uint8)
    
    if isinstance(graph, nx.DiGraph):
        graph = nx.Graph(graph)

    # Validate names and collect vocabulary
    node_names = {}
    for node, data in graph.nodes(data=True):
        if node_name_attr not in data:
            raise ValueError(f"Node {node!r} missing attribute {node_name_attr!r}!")
        name = str(data[node_name_attr])
        if name not in substitution_matrix:
            raise ValueError(f"Module name {name!r} not found in substitution matrix!")
        node_names[node] = name

    vocab = sorted(substitution_matrix.keys())
    vocab_index = {name: i for i, name in enumerate(vocab)}

    # Build substitution-aware module embeddings from the matrix rows, using truncated SVD
    S = np.zeros((len(vocab), len(vocab)), dtype=float)
    for a in vocab:
        for b in vocab:
            S[vocab_index[a], vocab_index[b]] = float(substitution_matrix[a].get(b, 0.0))

    # Column-center so very global offsets matter less
    S = S - S.mean(axis=0, keepdims=True)

    U, sigma, _ = np.linalg.svd(S, full_matrices=False)
    k = min(embedding_dim, len(sigma))
    if k == 0:
        raise ValueError("Could not derive embeddings from substitution matrix!")
    
    E = U[:, :k] * np.sqrt(np.maximum(sigma[:k], 0.0))
    embeddings = {name: E[vocab_index[name]].astype(np.float32) for name in vocab}
    
    # Start fingerprinting
    fp = np.zeros(n_bits, dtype=np.uint16 if counted else np.uint8)

    # Helpers
    def add_counted_feature(*parts: Any, weight: int = 1) -> None:
        idx = hash_to_index(*parts)
        fp[idx] += weight

    def add_binary_feature(*parts: Any, repeats: int = 1) -> None:
        for salt in range(repeats):
            idx = hash_to_index("salt", salt, *parts)
            fp[idx] = 1

    def add_feature(*parts: Any, weight: int = 1, repeats: int = 1) -> None:
        if counted:
            add_counted_feature(*parts, weight=weight)
        else:
            add_binary_feature(*parts, repeats=repeats)

    # initial node states = subsitution-aware embeddings
    states = {node: embeddings[name].copy() for node, name in node_names.items()}

    # Hash module-content features (bag of modules, substitution-aware)
    module_counts = Counter(node_names.values())
    for name, count in module_counts.items():
        vec_q = quantize(embeddings[name], decimals=2)

        # Hash exact count-aware content feature
        add_feature("module", name, count, vec_q, count, weight=1, repeats=1)

        # Hash a few count-independent bits too so multiplicity is not everything
        add_feature("module_presence", name, vec_q, weight=1, repeats=1)


    # # Hash edge-content features (helps compare adjacency while still being
    # # tolerant because names are represented by embeddings)
    # edge_counter = Counter()
    # for u, v in graph.edges():
    #     a = node_names[u]
    #     b = node_names[v]
    #     pair = tuple(sorted((a, b)))
    #     edge_counter[pair] += 1
    
    # for (a, b), count in edge_counter.items():
    #     va = embeddings[a]
    #     vb = embeddings[b]
    #     pair_vec = (va + vb) / 2.0
    #     pair_q = quantize(pair_vec, decimals=2)

    #     if count_similar_edges:
    #         add_feature("edge", a, b, count, pair_q, weight=1, repeats=1)

    #     add_feature("edge_presence", a, b, pair_q, weight=1, repeats=1)

    # # WL-like neighborhood propagation: captures topology, cycles, branching,
    # # disconnected modules, etc
    # for r in range(radius + 1):
    #     for node in graph.nodes():
    #         h = states[node]
    #         deg = graph.degree(node)
    #         name = node_names[node]

    #         # Hash local state
    #         add_feature("node_state", "r", r, "name", name, "deg", deg, "state", quantize(h, decimals=2), weight=1, repeats=1)

    #         # Also hash neighborhood composition a bit more explicitly
    #         nbr_names = sorted(node_names[nbr] for nbr in graph.neighbors(node))
    #         add_feature("node_nbr_names", "r", r, "name", name, "deg", deg, tuple(nbr_names), weight=1, repeats=1)

    #     if r == radius:
    #         break

    #     new_states = {}
    #     for node in graph.nodes():
    #         h_self = states[node]
    #         nbrs = list(graph.neighbors(node))

    #         if nbrs:
    #             nbr_vecs = np.stack([states[n] for n in nbrs], axis=0)
    #             h_nbr_mean = nbr_vecs.mean(axis=0)
    #             h_nbr_max = nbr_vecs.max(axis=0)
    #             deg_scale = 1.0 / math.sqrt(len(nbrs))
    #         else:
    #             h_nbr_mean = np.zeros_like(h_self)
    #             h_nbr_max = np.zeros_like(h_self)
    #             deg_scale = 1.0

    #         # Fixed deterministic update, not learned
    #         new_h = np.concatenate([
    #             h_self,
    #             deg_scale * h_nbr_mean,
    #             0.5 * h_nbr_max,
    #         ])

    #         # Compress back to embedding_dim deterministically by hashing chunks
    #         compressed = np.zeros(k, dtype=np.float32)
    #         for i, value in enumerate(new_h):
    #             j = i % k
    #             compressed[j] += float(value)

    #         # Nonlinearity + normalization for stability
    #         compressed = np.tanh(compressed)
    #         norm = np.linalg.norm(compressed)
    #         if norm > 0:
    #             compressed /= norm
            
    #         new_states[node] = compressed

    #     states = new_states

    # # Hash a few global graph features lightly
    # n_components = nx.number_connected_components(graph)
    # component_sizes = sorted(len(c) for c in nx.connected_components(graph))
    # n_nodes = graph.number_of_nodes()
    # n_edges = graph.number_of_edges()
    # cyclomatic = n_edges - n_nodes + n_components

    # add_feature("global", "components", n_components, tuple(component_sizes), weight=1, repeats=1)
    # add_feature("global", "nodes", n_nodes, "edges", n_edges, "cyclomatic", cyclomatic, weight=1, repeats=1)

    return fp
