class Graph {
    constructor() {
        this.adjList = {}; // { nodeId: [{ id: neighborId, bond: 'single' }, ...] }
    }

    /**
     * Adds a node to the graph.
     * @param {string} id - Unique identifier for the atom.
     */
    addNode(id) {
        if (!this.adjList[id]) {
            this.adjList[id] = [];
        }
    }

    /**
     * Adds an undirected edge (bond) between two nodes with a specified bond type.
     * @param {string} source - Source atom ID.
     * @param {string} target - Target atom ID.
     * @param {string} bond - Bond type (e.g., 'single', 'double', 'triple').
     */
    addEdge(source, target, bond = 'single') {
        this.addNode(source);
        this.addNode(target);
        this.adjList[source].push({ id: target, bond: bond });
        this.adjList[target].push({ id: source, bond: bond }); // Undirected graph
    }
}

class VF2Matcher {
    /**
     * Initializes the VF2Matcher with two graphs.
     * @param {Graph} graph1 - The first molecular graph.
     * @param {Graph} graph2 - The second molecular graph.
     */
    constructor(graph1, graph2) {
        this.graph1 = graph1;
        this.graph2 = graph2;
        this.mapping = {}; // Mapping from graph1 nodes to graph2 nodes
        this.mappedNodes1 = new Set();
        this.mappedNodes2 = new Set();
        this.result = null;
    }

    /**
     * Initiates the isomorphism check.
     * @returns {boolean} - True if isomorphic, else false.
     */
    isIsomorphic() {
        // Preliminary checks
        if (Object.keys(this.graph1.adjList).length !== Object.keys(this.graph2.adjList).length) {
            return false;
        }

        // Check degree sequence
        if (!this.checkDegreeSequence()) {
            return false;
        }

        return this.match(0);
    }

    /**
     * Checks if both graphs have identical degree sequences.
     * @returns {boolean}
     */
    checkDegreeSequence() {
        const degrees1 = Object.values(this.graph1.adjList)
            .map(neighbors => neighbors.length)
            .sort((a, b) => b - a);
        const degrees2 = Object.values(this.graph2.adjList)
            .map(neighbors => neighbors.length)
            .sort((a, b) => b - a);

        if (degrees1.length !== degrees2.length) return false;
        for (let i = 0; i < degrees1.length; i++) {
            if (degrees1[i] !== degrees2[i]) return false;
        }
        return true;
    }

    /**
     * Recursive matching function.
     * @param {number} depth - Current depth in the recursion.
     * @returns {boolean} - True if mapping is found, else false.
     */
    match(depth) {
        if (depth === Object.keys(this.graph1.adjList).length) {
            // All nodes are mapped
            this.result = { ...this.mapping };
            return true;
        }

        // Select the next node from graph1 to map
        const node1 = this.selectNode1(depth);
        const candidates = this.findCandidates(node1);

        for (let node2 of candidates) {
            if (this.isFeasible(node1, node2)) {
                // Extend the mapping
                this.mapping[node1] = node2;
                this.mappedNodes1.add(node1);
                this.mappedNodes2.add(node2);

                // Recurse
                if (this.match(depth + 1)) {
                    return true;
                }

                // Backtrack
                delete this.mapping[node1];
                this.mappedNodes1.delete(node1);
                this.mappedNodes2.delete(node2);
            }
        }

        return false; // No valid mapping found at this depth
    }

    /**
     * Selects the next node from graph1 to map based on a heuristic.
     * @param {number} depth - Current depth in the recursion.
     * @returns {string} - The selected node ID from graph1.
     */
    selectNode1(depth) {
        // Heuristic: select nodes in descending order of degree
        const nodes = Object.keys(this.graph1.adjList);
        nodes.sort((a, b) => this.graph1.adjList[b].length - this.graph1.adjList[a].length);
        return nodes[depth];
    }

    /**
     * Finds candidate nodes in graph2 that can be mapped to a given node in graph1.
     * @param {string} node1 - The node ID from graph1.
     * @returns {string[]} - Array of candidate node IDs from graph2.
     */
    findCandidates(node1) {
        const degree = this.graph1.adjList[node1].length;
        return Object.keys(this.graph2.adjList).filter(node2 =>
            this.graph2.adjList[node2].length === degree &&
            !this.mappedNodes2.has(node2)
        );
    }

    /**
     * Checks if mapping node1 to node2 is feasible.
     * @param {string} node1 - Node ID from graph1.
     * @param {string} node2 - Node ID from graph2.
     * @returns {boolean} - True if feasible, else false.
     */
    isFeasible(node1, node2) {
        const neighbors1 = this.graph1.adjList[node1];
        const neighbors2 = this.graph2.adjList[node2];

        for (let neighbor1 of neighbors1) {
            const mappedNeighbor2 = this.mapping[neighbor1.id];
            if (mappedNeighbor2) {
                // Find corresponding bond in graph2
                const bondInGraph2 = neighbors2.find(edge => edge.id === mappedNeighbor2);
                if (!bondInGraph2 || bondInGraph2.bond !== neighbor1.bond) {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * Retrieves the successful mapping after isomorphism is confirmed.
     * @returns {Object|null} - Mapping object or null if no mapping exists.
     */
    getMapping() {
        return this.result;
    }
}

export { Graph, VF2Matcher };