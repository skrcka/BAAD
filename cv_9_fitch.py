leaf_states = {
    'human': 'C',
    'chimp': 'T',
    'gibbon': 'G',
    'lemur': 'T',
    'gorilla': 'A',
    'bonobo': 'A',
}


class Node:
    def __init__(self, name: str) -> None:
        self.name = name
        self.children: list[Node] = []
        self.state: set[int] = set()
        self.assigned_state: dict[str, str] | None = None

    def add_child(self, child: 'Node') -> None:
        self.children.append(child)


def print_post_order(node: Node, level: int = 0) -> None:
    for child in node.children:
        print_post_order(child, level + 1)
    print('  ' * level + f"{node.name}: {node.state}")


def print_pre_order(node: Node, level: int = 0) -> None:
    print('  ' * level + f"{node.name}: {node.assigned_state}")
    for child in node.children:
        print_pre_order(child, level + 1)


def fitch_post_order(node: Node) -> set[str]:
    if not node.children:
        return node.state
    for child in node.children:
        node.state |= fitch_post_order(child)
    if len(node.children) == 2:
        intersection = node.children[0].state & node.children[1].state
        node.state = intersection if intersection else node.children[0].state | node.children[1].state
    return node.state


def fitch_pre_order(node: Node, parent_state: dict[str, str] | None = None) -> None:
    if parent_state and parent_state in node.state:
        node.assigned_state = parent_state
    elif node.state:
        node.assigned_state = next(iter(node.state))
    for child in node.children:
        fitch_pre_order(child, node.assigned_state)


def calculate_parsimony(node: Node, assigned_states: dict[str, str], parsimony_length: int = 0) -> int:
    for child in node.children:
        if child.assigned_state != node.assigned_state:
            parsimony_length += 1
        parsimony_length = calculate_parsimony(child, assigned_states, parsimony_length)
    if node.name not in leaf_states:
        assigned_states[node.name] = node.assigned_state
    return parsimony_length


def main() -> None:
    nodes = {
        'human': Node('human'),
        'chimp': Node('chimp'),
        'gibbon': Node('gibbon'),
        'lemur': Node('lemur'),
        'gorilla': Node('gorilla'),
        'bonobo': Node('bonobo'),
        'internode1': Node('internode1'),
        'internode2': Node('internode2'),
        'internode3': Node('internode3'),
        'internode4': Node('internode4'),
        'root': Node('root')
    }

    nodes['internode1'].add_child(nodes['human'])
    nodes['internode1'].add_child(nodes['chimp'])
    nodes['internode2'].add_child(nodes['gibbon'])
    nodes['internode2'].add_child(nodes['lemur'])
    nodes['internode3'].add_child(nodes['internode2'])
    nodes['internode3'].add_child(nodes['gorilla'])
    nodes['internode4'].add_child(nodes['internode1'])
    nodes['internode4'].add_child(nodes['internode3'])
    nodes['root'].add_child(nodes['internode4'])
    nodes['root'].add_child(nodes['bonobo'])

    for leaf, state in leaf_states.items():
        nodes[leaf].state.add(state)

    fitch_post_order(nodes['root'])
    fitch_pre_order(nodes['root'])
    assigned_states = {}
    parsimony_length = calculate_parsimony(nodes['root'], assigned_states)

    print("Uzly po prechode zdole-hore:")
    print_post_order(nodes['root'])
    print("\nUzly po prechode zhora-dole")
    print_pre_order(nodes['root'])
    print("\nParsimonie:", parsimony_length)
    print("\nStavy vnutornych uzlov", assigned_states)


if __name__ == '__main__':
    main()