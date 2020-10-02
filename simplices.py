# simplices.py
# store simplicial complexes and compute homology
# (c) 2020, Samuel Rabinowitz

from itertools import chain, combinations

# gives a where len(a_set) = 2^a

def dim_set(a_set):
    return len(a_set).bit_length() - 1

# individual simplex, immutable set of vertices

class Simplex(frozenset):

    def __add__(self, other):
        return SimplexSum([self]) + other

    def __radd__(self, other):
        return self.__add__(other)

    def __str__(self):
        if self.is_zero():
            return "∅"
        else:
            return "".join([str(v) for v in self])

    def order(self):
        return len(self) - 1

    # checks if empty simplex

    def is_zero(self):
        return (self.order() < 0)

    # returns all subsets of self as simplices
    # save the empty simplex

    def get_sub_simplices(self):
        l_v = list(self)
        order = self.order()
        for verts in chain.from_iterable(
            combinations(l_v, r) for r in
            range(1, order + 1)):
            yield((Simplex(verts), len(verts) - 1))
        yield((self, order))

    # simplex boundary

    def boundary(self):
        l_v = list(self)
        order = self.order()
        if(order < 1):
            return EmptySimplexSum
        else:
            return SimplexSum([Simplex(verts) for verts
                in combinations(l_v, order)])

# sum of simplices; user can should add regular simplices
# together to get a SimplexSum instead of trying to initialize
# one directly

class SimplexSum(frozenset):

    def __add__(self, other):
        if(other == 0):
            return self
        if(isinstance(other, Simplex)):
            to_add = SimplexSum([other])
        elif(isinstance(other, SimplexSum)):
            to_add = other
        else:
            raise TypeError(
                "Expected Simplex or SimplexSum, got {}.".format(
                    type(other)))
        return SimplexSum((self ^ to_add) - EmptySimplexSum)

    def __radd__(self, other):
        return self.__add__(other)

    def __str__(self):
        return " + ".join([str(simplex) for simplex in self])

    def boundary(self):
        bd = sum([simplex.boundary() for simplex in self])
        if(bd == 0 or len(bd) == 0):
            bd = EmptySimplexSum
        return bd

    # checks if empty; len(self) technically should never be 0
    # because users should always do EmptySimplexSum to get an empty
    # simplex sum instead of trying to initializing one themselves

    def is_zero(self):
        return (len(self) == 0 or (len(self) == 1
            and list(self)[0].is_zero()))

# allows for easy printing of set of Simplex's or
# SimplexSum's

class SimplexList(set):

    def __str__(self):
        return "{" + ", ".join([str(v) for v in self]) + "}"

# mutable class for storing simplices,
# performing calculations on data

class SimplicialComplex:

    def __init__(self):
        self.__simplices = []

    def __repr__(self):
        return "SimplicialComplex({})".format(self.__simplices)

    def __str__(self):
        return "{" + ", ".join([str(simplex) for order in
            self.__simplices for simplex in order]) + "}"

    # add new simplex
    # note: if we add a higher order simplex,
    # all lower order sub-simplices are added
    # e.g. adding ABC also adds A, B, C, AB,
    # AC, BC

    def add_simplex(self, simplex):
        simplex = Simplex(simplex)

        k = simplex.order()
        mx_s = len(self.__simplices)
        for _ in range(mx_s, k + 1):
            self.__simplices.append(SimplexList())

        for s_simplex, s_k in simplex.get_sub_simplices():
            self.__simplices[s_k].add(s_simplex)

    # provide iterable list of simplices

    def add_simplices(self, simplices):
        for simplex in simplices:
            self.add_simplex(simplex)

    # get all k-simplices

    def get_k_simplices(self, k = 0):
        if(k < 0 or k >= len(self.__simplices)):
            return SimplexList()
        else:
            return self.__simplices[k]

    # chain group of order k

    def C_k(self, k):
        if(k < 0 or k >= len(self.__simplices)):
            return SimplexList([EmptySimplexSum])
        l_k_s = list(self.get_k_simplices(k = k))
        l = len(l_k_s)
        return SimplexList([SimplexSum(simplices) for simplices in
            chain.from_iterable(combinations(l_k_s, r) for r in
                range(1, l + 1))] + [EmptySimplexSum])

    # cycle group of order k

    def Z_k(self, k):
        return SimplexList([z for z in self.C_k(k)
            if z.boundary().is_zero()])

    # boundary group of order k

    def B_k(self, k):
        return SimplexList([c.boundary() for c in self.C_k(k + 1)])

    # homology group of order k

    def H_k(self, k):
        z_k = self.Z_k(k)
        b_k = self.B_k(k)
        h_k = SimplexList()
        for ss_test in z_k:
            if not any([((ss_test + ss_base) in b_k)
                for ss_base in h_k]):
                h_k.add(ss_test)
        return h_k

    # Betti number of order k

    def beta_k(self, k):
        return dim_set(self.H_k(k))

# default empty simplices;
# use these, do not define your own

EmptySimplex = Simplex()
EmptySimplexSum = SimplexSum([EmptySimplex])