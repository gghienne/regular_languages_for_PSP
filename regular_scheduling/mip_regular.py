from docplex.mp.model import Model
from regular_scheduling.operations import intersection,unfold
from collections import defaultdict
from itertools import product

from regular_scheduling.tools import show

class Regular_Model(Model):

    def __init__(self,name = None,context = None,**kwargs):
        super().__init__(name,context,**kwargs)
        self.regular_solutions = {}
        self.building_time = None

    def add_regular_constraint(self,x,number,sequence,*dfas, ctname = None):

        dfa = intersection(*dfas)
        dfa = unfold(dfa,sequence)

        if type(sequence) == int:
            sequence = [dfa["alphabet"]]*sequence
        T = len(sequence)

        delta_minus,delta_plus,sigma = [defaultdict(set) for _ in range(T)],[defaultdict(set) for _ in range(T)],[defaultdict(set) for _ in range(T)]

        level = {dfa["initial_state"]}
        for t in range(T):
            for q1,s in [(q,s) for q,s in product(level,sequence[t]) if (q,s) in dfa["transitions"]]:
                q2 = dfa["transitions"][q1,s]
                delta_plus[t][q1].add((s,q2))
                delta_minus[t][q2].add((q1,s))
                sigma[t][s].add((q1,q2))
            level = delta_minus[t].keys()

        start = min(t for t,_ in x.keys())

        f = self.continuous_var_dict(((q1,s,q2) for c in delta_plus for q1,arcs in c.items() for s,q2 in arcs),lb = 0, ub = number)   
        
        flow_constraints = []
        flow_constraints += [self.add_constraint(self.sum(f[q0,s,q] for q0,arcs in delta_plus[0].items() for s,q in arcs) == number)]        
        flow_constraints += self.add_constraints(self.sum(f[q1,s,q] for q1,s in delta_minus[i][q]) == self.sum(f[q,s,q2] for s,q2 in c[q]) for i,c in enumerate(delta_plus[1:]) for q in c)
        flow_constraints += self.add_constraints(x[t,s] == self.sum(f[q1,s,q2] for q1,q2 in sigma[t-start][s]) for t,s in x.keys())

        def post_processing():
            if not self.solution:
                return None
            regular_solution = [[None]*T for _ in range(number)]
            f_star = {(q1,s,q2):round(self.solution.get_value(f[q1,s,q2])) for q1,s,q2 in f.keys()}
            for i in range(number):   
                q1 = dfa["initial_state"]
                for t in range(T):
                    s,q2 = next((s,q) for s,q in delta_plus[t][q1] if f_star[q1,s,q])
                    f_star[q1,s,q2] -= 1
                    regular_solution[i][t],q1 = s,q2
            return regular_solution
        
        if ctname:
            self.regular_solutions[ctname] = post_processing
        self.regular_solutions[str(flow_constraints)] = post_processing

        return flow_constraints
    
    def add_regular_constraint_(self,x,number,sequence,*dfas, ctname = None):

        self.add_regular_constraint(self,x,number,sequence,*dfas, ctname = ctname)
    
    def get_regular_solution(self,regular_ct):

        return self.regular_solutions[str(regular_ct)]()
