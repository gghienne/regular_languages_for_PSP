from PySimpleAutomata import DFA
from collections import defaultdict
from itertools import product

def minimize(dfa):
    """Returns the minimal DFA with the same language as input DFA."""
    return _rename(DFA.dfa_trimming(DFA.dfa_minimization(_rename(dfa))))


def intersection(*dfas):
    """Returns the intersection of the DFAs in input."""
    dfa1 = dfas[0]
    for dfa2 in dfas[1:]:
        dfa1 = DFA.dfa_intersection(dfa1, dfa2)
    alphabet = dfa1["alphabet"]
    dfa = _rename(DFA.dfa_trimming(dfa1))
    dfa["alphabet"] = alphabet
    return dfa


def union(*dfas):
    """Returns the union of the DFAs in input."""
    dfa1 = dfas[0]
    for dfa2 in dfas[1:]:
        dfa1 = DFA.dfa_union(dfa1, dfa2)
    alphabet = dfa1["alphabet"]
    dfa = _rename(DFA.dfa_trimming(dfa1))
    dfa["alphabet"] = alphabet
    return dfa


def complementation(dfa):
    """Returns the complement of the DFA in input."""
    new_dfa = _rename(DFA.dfa_trimming(DFA.dfa_complementation(dfa)))
    new_dfa["alphabet"] = dfa["alphabet"]
    return new_dfa


def side(dfa, left = "", right = ""):
    """Returns the result of the corresponding side-hypothesis operation."""
    dfa_copy = dfa.copy()
    dfa_copy["initial_state"] = __run(left, dfa_copy)
    finals = set()
    for q in dfa_copy["states"]:
        if __run(right, dfa_copy, q) in dfa_copy["accepting_states"]:
            finals.add(q)
    dfa_copy["accepting_states"] = finals
    return dfa_copy


def windows(dfa, period, lower_bound = 0, upper_bound = float("inf")):
    """Returns the result of the corresponding windows-cardinality operation."""

    dfa = unfold(dfa, period)

    alphabet_dict = {}
    for s in dfa["alphabet"]:
        key = tuple(dfa["transitions"].get((q, s)) for q in dfa["states"])
        alphabet_dict.setdefault(key, tuple())
        alphabet_dict[key] += (s,)
    alphabet_dict = {v[0]: v for v in alphabet_dict.values()}

    new_dfa = __init_dfa(alphabet_dict.keys(), (tuple(), 0))

    OPEN = {new_dfa["initial_state"]}
    while OPEN:
        state = OPEN.pop()
        new_dfa["states"].add(state)
        w, n = state
        if lower_bound <= n:
            new_dfa["accepting_states"].add(state)
        if upper_bound == float("inf") and n == lower_bound:
            new_dfa["transitions"].update(
                {
                    (state, t): state
                    for s in new_dfa["alphabet"]
                    for t in alphabet_dict[s]
                }
            )
        else:
            for s in new_dfa["alphabet"]:
                v = (w + (s,))[1 - period :]
                while __run(v, dfa) is None:
                    v = v[1:]
                if (
                    len(w) == period - 1
                    and __run(w + (s,), dfa) in dfa["accepting_states"]
                ):
                    m = n + 1
                else:
                    m = n
                if m <= upper_bound:
                    for t in alphabet_dict[s]:
                        if upper_bound == float("inf") and m == lower_bound:
                            new_dfa["transitions"][state, t] = (tuple(), m)
                        else:
                            new_dfa["transitions"][state, t] = (v, m)
                    if new_dfa["transitions"][state, s] not in new_dfa["states"]:
                        OPEN.add(new_dfa["transitions"][state, s])

    new_dfa["alphabet"] = set().union(*alphabet_dict.values())
    return _rename(new_dfa)


def periodic(dfa, period, binary_dfa = None):
    """Returns the result of the corresponding periodic-composition operation."""

    if not binary_dfa:
        binary_dfa = ones_dfa
    dfa = unfold(dfa, period)
    new_dfa = __init_dfa(
        dfa["alphabet"], (dfa["initial_state"], binary_dfa["initial_state"], 0)
    )

    OPEN = {new_dfa["initial_state"]}
    while OPEN:
        state = OPEN.pop()
        new_dfa["states"].add(state)
        q, q_01, l = state
        if q_01 in binary_dfa["accepting_states"]:
            new_dfa["accepting_states"].add(state)
        for s in dfa["alphabet"]:
            if l == period - 1:
                s_01 = (
                    "1"
                    if dfa["transitions"].get((q, s)) in dfa["accepting_states"]
                    else "0"
                )
                next_q_01 = binary_dfa["transitions"].get((q_01, s_01))
                if next_q_01 is not None:
                    new_dfa["transitions"][state, s] = (
                        dfa["initial_state"],
                        next_q_01,
                        0,
                    )
                    if new_dfa["transitions"][state, s] not in new_dfa["states"]:
                        OPEN.add(new_dfa["transitions"][state, s])
            else:
                next_q = dfa["transitions"].get((q, s))
                new_dfa["transitions"][state, s] = (next_q, q_01, l + 1)
                if new_dfa["transitions"][state, s] not in new_dfa["states"]:
                    OPEN.add(new_dfa["transitions"][state, s])

    return _rename(new_dfa)


def mask(dfa, binary_mask):
    """Returns the result of the corresponding periodic-mask operation."""
    new_dfa = __init_dfa(dfa["alphabet"], (dfa["initial_state"], 0))

    OPEN = {new_dfa["initial_state"]}
    while OPEN:
        state = OPEN.pop()
        new_dfa["states"].add(state)
        q, l = state
        if q in dfa["accepting_states"]:
            new_dfa["accepting_states"].add(state)
        for s in dfa["alphabet"]:
            if binary_mask[l] == "1":
                next_q = dfa["transitions"].get((q, s))
            else:
                next_q = q
            if next_q is not None:
                new_dfa["transitions"][state, s] = (next_q, (l + 1) % len(binary_mask))
                if new_dfa["transitions"][state, s] not in new_dfa["states"]:
                    OPEN.add(new_dfa["transitions"][state, s])
    return _rename(new_dfa)

def unfold(dfa, sequence):
    """Returns the corresponding unfolded DFA."""
    if type(sequence) == int:
        sequence = [dfa["alphabet"]]*sequence
    
    # initialize
    n = len(sequence)
    current_layer = {dfa["initial_state"]}
    children = [defaultdict(set) for _ in range(n)]

    # expand forward
    for i in range(n):
        next_layer = set()
        for q1,s in [(q,s) for q,s in product(current_layer,sequence[i]) if (q,s) in dfa["transitions"]]:
            q2 = dfa["transitions"][q1,s]
            children[i][q1].add((s,q2))
            next_layer.add(q2)
        current_layer = next_layer

    # filter and merge last layer
    final_states = tuple(current_layer.intersection(dfa["accepting_states"]))
    new_state = {q:final_states for q in final_states}
    if not new_state:
        return  {
            "alphabet": dfa["alphabet"],
            "states": {1},
            "initial_state": 1,
            "transitions": {},
            "accepting_states": {},
        }
    
    # validate backward
    delta = {}
    for i in range(n-1,-1,-1):
        for q1,arcs in children[i].items():
            children[i][q1] = tuple(sorted((s,new_state[q2]) for s,q2 in arcs if q2 in new_state))
        reverse_children,new_state = defaultdict(list),{}
        for q,arcs in children[i].items():
            reverse_children[arcs].append(q)
        for arcs,states in reverse_children.items():
            states = tuple(states)
            if arcs:
                children[i][states] = arcs
                for s,q in arcs:
                    delta[(i,states),s] = (i+1,q)
                for q in states:
                    new_state[q] = states
            for q in states:
                del children[i][q]

    unfolded = {
            "alphabet": dfa["alphabet"],
            "states": {(i,q) for i,c in enumerate(children) for q in c}.union({(n,final_states)}),
            "initial_state": (0,(dfa["initial_state"],)),
            "transitions": delta,
            "accepting_states": {(n,final_states)},
        }
    return _rename(unfolded)

def __init_dfa(alphabet, initial_state):
    return {
        "alphabet": alphabet,
        "states": set(),
        "initial_state": initial_state,
        "transitions": {},
        "accepting_states": set(),
    }


def __run(word, dfa, state=None):
    if state is None:
        state = dfa["initial_state"]
    for s in word:
        state = dfa["transitions"].get((state, s))
        if state is None:
            return None
    return state


def _rename(dfa):
    Q = {dfa["initial_state"]: 1}
    Q.update(
        {
            q: i + 2
            for i, q in enumerate(
                dfa["states"].difference(
                    dfa["accepting_states"].union({dfa["initial_state"]})
                )
            )
        }
    )
    Q.update(
        {
            q: len(dfa["states"]) - i
            for i, q in enumerate(
                dfa["accepting_states"].difference({dfa["initial_state"]})
            )
        }
    )
    return {
        "alphabet": dfa["alphabet"],
        "states": set(range(1, len(Q) + 1)),
        "initial_state": Q[dfa["initial_state"]],
        "transitions": {(Q[k[0]], k[1]): Q[v] for k, v in dfa["transitions"].items()},
        "accepting_states": set(Q[q] for q in dfa["accepting_states"]),
    }

ones_dfa = {
        "alphabet": {"0","1"},
        "states": {0},
        "initial_state": 0,
        "transitions": {(0,"1"):0},
        "accepting_states": {0},
    }