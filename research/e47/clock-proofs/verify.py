"""Exact block spectra, kernel ranks and OP-4 endpoints. Run with SymPy."""
from itertools import product
from collections import defaultdict
import sympy as sp

fibers = defaultdict(list)
for triple in product(range(5), repeat=3):
    fibers[sum(triple) - 6].append(triple)

x = sp.symbols("x")
multiplicities = [1, 3, 5, 4, 3, 2, 1]
nullities = []
for weight, basis in sorted(fibers.items()):
    identity = sp.eye(len(basis))
    casimir = weight * (weight + 1) * identity
    for column, state in enumerate(basis):
        for i in range(3):
            if state[i] == 4:
                continue
            raised = list(state)
            raised[i] += 1
            for j in range(3):
                if raised[j] > 0:
                    target = raised.copy()
                    target[j] -= 1
                    casimir[basis.index(tuple(target)), column] += (4 - state[i]) * raised[j]
    expected = sp.prod((x-j*(j+1))**multiplicities[j] for j in range(abs(weight), 7))
    assert casimir.charpoly(x).as_poly() == sp.Poly(expected, x)
    selector = (casimir - 6*identity) * (casimir - 30*identity)
    nullities.append(len(basis) - selector.rank())
assert nullities == [0, 2, 2, 2, 7, 7, 7, 7, 7, 2, 2, 2, 0]
assert sum(nullities) == 47

t = sp.symbols("t", nonnegative=True)
curve = 1 - sum(sp.Rational(d,125)*sp.exp(-2*(j*(j+1))**2*t)
                for j,d in [(0,1),(1,9),(3,28),(4,27),(6,13)])
assert curve.subs(t,0) == sp.Rational(47,125)
assert sp.limit(curve,t,sp.oo) == sp.Rational(124,125)
assert curve.diff(t).is_positive
print("PASS: 13 exact spectra, kernel dimension 47, endpoints and positive derivative")
