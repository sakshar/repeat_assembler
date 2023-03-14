cycle_start, cycle_end = -1, -1


def dfs(v, color, parent, adj):
    global cycle_start, cycle_end
    color[v] = 1
    for u in adj[v]:
        if color[u] == 0:
            parent[u] = v
            if dfs(u, color, parent, adj):
                return True
        elif color[u] == 1:
            cycle_end = v
            cycle_start = u
            return True
    color[v] = 2
    return False


def find_cycle(n, adj):
    nodes = list(adj.keys())
    color = dict()
    parent = dict()
    for node in nodes:
        color[node] = 0
        parent[node] = -1
    global cycle_start, cycle_end
    cycle_start = -1

    for v in nodes:
        if color[v] == 0 and dfs(v, color, parent, adj):
            break
    cycle = list()
    if cycle_start == -1:
        print("Acyclic")
    else:
        cycle.append(cycle_start)
        v = cycle_end
        while v != cycle_start:
            cycle = [v] + cycle
            v = parent[v]
        cycle = [cycle_start] + cycle

        print("Cycle found: ", cycle)
    return cycle


n = 4
adj = dict()
adj[0] = [1, 2]
adj[1] = [2]
adj[2] = [4]
adj[4] = [0]
print(find_cycle(n, adj))
