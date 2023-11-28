distance_matrix = {
    'A': {'A': 0, 'B': 20, 'C': 26, 'D': 26, 'E': 26},
    'B': {'A': 20, 'B': 0, 'C': 26, 'D': 26, 'E': 26},
    'C': {'A': 26, 'B': 26, 'C': 0, 'D': 16, 'E': 16},
    'D': {'A': 26, 'B': 26, 'C': 16, 'D': 0, 'E': 10},
    'E': {'A': 26, 'B': 26, 'C': 16, 'D': 10, 'E': 0},
}


def compute_new_distances(matrix: dict[str, dict[str, int]], a: str, b: str) -> dict[str, int]:
    return {node: (matrix[a][node] + matrix[b][node]) / 2 for node in matrix if node != a and node != b}


def main() -> None:
    leaf_nodes = list(distance_matrix.keys())
    node_distances = {}

    while len(leaf_nodes) > 1:
        min_distance = float('inf')
        min_i = None
        min_j = None

        for i, leaf_node1 in enumerate(leaf_nodes):
            for leaf_node2 in leaf_nodes[i+1:]:
                dist = distance_matrix[leaf_node1][leaf_node2]
                if dist < min_distance:
                    min_distance = dist
                    min_i, min_j = leaf_node1, leaf_node2

        new_node_label = f"({min_i}, {min_j})"
        new_distances = compute_new_distances(distance_matrix, min_i, min_j)

        distance_matrix[new_node_label] = new_distances
        leaf_nodes.append(new_node_label)

        node_distances[new_node_label] = min_distance

        del distance_matrix[min_i]
        del distance_matrix[min_j]

        for node, distances in distance_matrix.items():
            if node != new_node_label:
                distances[new_node_label] = new_distances[node]

        leaf_nodes = list(filter(lambda node: node != min_i and node != min_j, leaf_nodes))

    phylogenetic_tree = leaf_nodes[0]
    print("Fylogenetic tree (UPGMA) with distances:")
    print(phylogenetic_tree)
    print(f"Distance from root: {node_distances[phylogenetic_tree]}")


if __name__ == '__main__':
    main()
