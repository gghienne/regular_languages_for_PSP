from bisect import bisect
from collections import defaultdict
from itertools import product
from regular_scheduling.operations import intersection,complementation,_rename

def cardinality(symbols,lower_bounds = tuple(),upper_bounds = tuple()):
    """ Return the corresponding cardinality rule."""
    symbols = [tuple(k) if type(k)!=str else (k,) for k in symbols]
    if not lower_bounds:
        lower_bounds = tuple([0]*len(symbols))
    if not upper_bounds:
        upper_bounds = tuple([float("inf")]*len(symbols))
    lower = {k:lower_bounds[i] for i,k in enumerate(symbols)}
    upper = {k:upper_bounds[i] for i,k in enumerate(symbols)}  
    maxi = {k:v if v<float("inf") else lower[k] for k,v in upper.items()}
    states = set(product(*[range(v+1) for v in maxi.values()]))
    transitions = {(q,s):q[:i] + (int(q[i]+int(q[i] < maxi[k])),) + q[i+1:]
                  for q in states for i,(k,v) in enumerate(upper.items()) for s in k if q[i]<v}
    finals = set(product(*[range(v,maxi[k]+1) for k,v in lower.items()]))
    dfa = {"alphabet"         :set(s for k in symbols for s in k),
           "states"           :states,
           "initial_state"    :tuple([0]*len(lower)),
           "transitions"      :transitions,
           "accepting_states" :finals}
    
    return _rename(dfa)

def stretch(symbols,lower_bounds = tuple(),upper_bounds = tuple()):
    """ Return the corresponding stretch rule."""
    symbols = [tuple(k) if type(k)!=str else (k,) for k in symbols]
    if not lower_bounds:
        lower_bounds = tuple([1]*len(symbols))
    if not upper_bounds:
        upper_bounds = tuple([float("inf")]*len(symbols))
    bounded_idx,unbounded_idx = [],[]
    for i,l in enumerate(lower_bounds):
        if (l,upper_bounds[i]) != (1,float("inf")):
            bounded_idx.append(i)
        else:
            unbounded_idx.append(i)
    lower = {symbols[i]:lower_bounds[i] for i in bounded_idx}
    upper = {symbols[i]:upper_bounds[i] for i in bounded_idx} 
    maxi = {k:v if v<float("inf") else lower[k] for k,v in upper.items()}
    states = set((k,m) for k,v in maxi.items() for m in range(1,v+1)).union({(("",),0)})
    finals = set((k,m) for k,v in maxi.items() for m in range(lower[k],v+1)).union({(("",),0)})
    transitions = {**{((k1,m),s):(k2,m*(k1==k2)+int(m*(k1==k2)<v))
                      for k1,m in states for k2,v in maxi.items() for s in k2
                      if ((m*(k1==k2) < upper[k2]) and (((k1,m) in finals) or (k1==k2)))},
                   **{((("",),0),s):(("",),0) for i in unbounded_idx for s in symbols[i]},
                   **{(q,s):(("",),0) for i in unbounded_idx for s in symbols[i] for q in finals}}
    dfa = {"alphabet"         :set(s for k in symbols for s in k),
           "states"           :states,
           "initial_state"    :(("",),0),
           "transitions"      :transitions,
           "accepting_states" :finals}
    return _rename(dfa)

def pattern(symbols,sequence,lower_bounds = tuple(),upper_bounds = tuple()):
    """ Return the corresponding patterned-stretch rule."""
    symbols = [tuple(k) if type(k)!=str else (k,) for k in symbols]
    sequence = [tuple(k) if type(k)!=str else (k,) for k in sequence]
    if not lower_bounds:
        lower_bounds = tuple([1]*len(sequence))
    if not upper_bounds:
        upper_bounds = tuple([float("inf")]*len(sequence))
    maxi = [u if u<float("inf") else lower_bounds[i] for i,u in enumerate(upper_bounds)]
    states = set((i,u) for i,m in enumerate(maxi) for u in range(1,m+1)).union({(-1,float("inf"))})
    finals = set((len(maxi)-1,m) for m in range(lower_bounds[len(maxi)-1],maxi[-1]+1))
    transitions = {**{((i,m),s):(i+1,1) for i,m in states if i<len(maxi)-1 and m>=lower_bounds[i] for s in tuple(sequence[i+1])},
                   **{((i,m),s):(i,m+int(m<maxi[i])) for i,m in states for s in tuple(sequence[i]) if m<upper_bounds[i]}}
    dfa = {"alphabet"         :set(s for k in symbols for s in k),
           "states"           :states,
           "initial_state"    :(-1,float("inf")),
           "transitions"      :transitions,
           "accepting_states" :finals}
    return _rename(dfa)

def knapsack(symbols,weights,lower_bound = 0,upper_bound = float("inf")):
    """ Return the corresponding knapsack rule."""
    symbols = [tuple(k) if type(k)!=str else (k,) for k in symbols]
    if lower_bound == 0:
        if upper_bound == float("inf"):
            return __universal(symbols)
        return __knapsack(symbols,weights,upper_bound)
    return intersection(
                    __knapsack(symbols,weights,upper_bound),
                    complementation(__knapsack(symbols,weights,lower_bound-1e-10)))

def __universal(symbols):
    symbols = [tuple(k) if type(k)!=str else (k,) for k in symbols]
    return {"alphabet"         :set(s for k in symbols for s in k),
            "states"           :{0},
            "initial_state"    :{0},
            "transitions"      :{(0,s):0 for k in symbols for s in k},
            "accepting_states" :{0}}

def __knapsack(symbols,weights,upper_bound):
    sorted_weights = sorted(set(weights).difference({0}))
    OPEN,CLOSED,MAX = [0],[],defaultdict(int)
    max_i = len(sorted_weights)
    while max_i > 0:
        q1 = OPEN.pop(0)
        CLOSED.append(q1)
        MAX[q1] = 1
        while MAX[q1] <= max_i:
            q2 = round(q1+sorted_weights[MAX[q1]-1],5)
            if q2 <= upper_bound:
                i = bisect(OPEN,q2)
                if i == 0 or OPEN[i-1] != q2:
                    OPEN.insert(i,q2)
                MAX[q1] += 1
            else:
                max_i = MAX[q1]-1
        MAX[q1] = max_i
    states,CLASS = set(),{}    
    last_merge,current_class = sorted_weights[0],CLOSED.pop()
    OPEN.append(current_class)
    CLASS = {q:current_class for q in OPEN}
    while round(current_class+sorted_weights[0],5) > last_merge:
        sigma = CLOSED.pop()
        i = 0
        if MAX[sigma] == MAX[current_class]:
            while i < MAX[sigma]:
                i += 1
                if CLASS[round(current_class+sorted_weights[i-1],5)] !=CLASS[round(sigma+sorted_weights[i-1],5)]:
                    i = MAX[sigma]+1
        if i == MAX[sigma]:
            last_merge = sigma
            OPEN.append(sigma)
        else:
            current_class = sigma
            states.add(tuple(OPEN))
            OPEN = [sigma]
        CLASS[sigma] = current_class
    states.update([tuple(OPEN)]+[(q,) for q in CLOSED])
    CLASS = {sigma:q for q in states for sigma in q}
    transitions = {(q,s):CLASS[round(q[0]+weights[i],5)] for q in states
                   for i,S in enumerate(symbols) for s in tuple(S)
                   if weights[i] <= ([0]+sorted_weights)[MAX[q[0]]]}
    return {"alphabet":         set(s for S in symbols for s in tuple(S)),
            "states":           states,
            "initial_state":    (0,),
            "transitions":      transitions,
            "accepting_states": states.copy()}