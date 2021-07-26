import copy

def check_G(G):
    for v in G:
        for v2 in G[v]:
            if v not in G[v2]:return False
    return True

def sort_G(G):
    return {v:sorted(G[v]) for v in sorted(G)}

def get_r(n, offset=0):
    G = {}
    for i in range(n):
        G[i + offset] = []
        G[i + offset].append((i + 1) % n + offset)
        G[i + offset].append((i - 1 + n) % n + offset)
    return G

def get_c(n, offset = 0):
    G = {offset:[offset+1]}
    for i in range(1, n-1):
        G[i + offset] = []
        G[i + offset].append(i + offset + 1)
        G[i + offset].append(i + offset -1)
    G[n-1 + offset] = [n-2 + offset]
    return G

def get_rc(Nr, Nc, Nrc):
    """
    generate ring with a chain 
    Nr - ring size
    Nc - chain size (outside of the ring)
    Nrc - number of atoms that overlap (at least 1)
    """
    G = comb_G(get_r(Nr), get_c(Nc + Nrc, -Nc))
    Gc2 = get_c(Nc + Nrc, Nr - Nrc)
    replace_v_in_G(Gc2, [(Nr - Nrc + i,i) for i in range(Nrc)])
    G = comb_G(G, Gc2)
    return sort_G(G)

def comb_G(G1, G2):
    G = copy.deepcopy(G1)
    for v in G2:
        if v not in G:
            G[v] = []
        for vv in G2[v]:
            if vv not in G[v]:
                G[v].append(vv)
    return G

def connect_v(G, v1, v2):
    if v2 not in G[v1]:
        G[v1].append(v2)
    if v1 not in G[v2]:
        G[v2].append(v1)

def replace_v_in_G(temp_G, v_pairs):
    for v_pair in v_pairs:
        temp_G[v_pair[1]] = temp_G.pop(v_pair[0])
        for v in temp_G:
            if v_pair[0] in temp_G[v]:
                temp_G[v].remove(v_pair[0])
                temp_G[v].append(v_pair[1])

