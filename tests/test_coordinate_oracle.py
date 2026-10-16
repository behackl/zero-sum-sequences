"""Cross-checks against independent coordinate arithmetic on small groups.

The oracles below sum integer tuples modulo the cyclic orders directly and
enumerate multisets with ``itertools``; they never use the package's own
sequence predicates as ground truth.
"""

from collections import Counter
from functools import cache
from itertools import combinations_with_replacement, product
from math import inf, lcm

import networkx as nx
import pytest

from zero_sum_sequences import (
    AdditiveSequenceSpace,
    AtomCatalogue,
    FiniteAdditiveGroup,
    zero_sum_constant,
)

GROUPS = [(1,), (2,), (3,), (4,), (5,), (6,), (2, 2), (2, 3), (2, 4), (2, 2, 2), (3, 3)]


def elements_of(moduli):
    return tuple(product(*(range(n) for n in moduli)))


def vector_sum(terms, moduli):
    return tuple(sum(t[i] for t in terms) % n for i, n in enumerate(moduli))


def sub_multisets(terms):
    counts = Counter(terms)
    keys = sorted(counts)
    for chosen in product(*(range(counts[k] + 1) for k in keys)):
        yield tuple(k for k, c in zip(keys, chosen) for _ in range(c))


def is_atom(terms, moduli):
    zero = (0,) * len(moduli)
    return bool(terms) and vector_sum(terms, moduli) == zero and not any(
        0 < len(sub) < len(terms) and vector_sum(sub, moduli) == zero
        for sub in sub_multisets(terms)
    )


def multisets(elements, max_length):
    for length in range(max_length + 1):
        yield from combinations_with_replacement(elements, length)


@pytest.mark.parametrize("moduli", GROUPS, ids=str)
def test_predicates_and_factorizations(moduli):
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(*moduli))
    elements = elements_of(moduli)
    zero = (0,) * len(moduli)
    nonzero = tuple(t for t in elements if t != zero)
    inputs = list(multisets(elements, 5 if len(elements) <= 8 else 4))
    atoms = [t for t in inputs if zero not in t and is_atom(t, moduli)]
    catalogue = AtomCatalogue(space, [space(a) for a in atoms])
    cache_ = space.factorization_cache(catalogue)
    vectors = [tuple(a.count(t) for t in nonzero) for a in atoms]
    atom_ids = {space(a): i for i, a in enumerate(atoms)}

    @cache
    def factorizations(state, start=0):
        """Nondecreasing tuples of atom indices whose vectors sum to ``state``."""

        if not any(state):
            return ((),)
        found = []
        for i in range(start, len(vectors)):
            if all(c >= d for c, d in zip(state, vectors[i])):
                rest = tuple(c - d for c, d in zip(state, vectors[i]))
                found.extend((i, *f) for f in factorizations(rest, i))
        return tuple(found)

    for terms in inputs:
        sequence = space(terms)
        assert tuple(sequence.total()) == vector_sum(terms, moduli)
        assert sequence.is_atom() == is_atom(terms, moduli)
        if zero in terms:
            continue

        expected = factorizations(tuple(terms.count(t) for t in nonzero))
        lengths = {len(f) for f in expected}
        local = sequence.factorization_solver()
        indexed = sequence.factorization_solver(atom_catalogue=catalogue)
        assert local.minimum_factorization_length() == min(lengths, default=None)
        assert local.maximum_factorization_length() == max(lengths, default=None)
        assert local.length_set() == indexed.length_set() == lengths
        assert cache_.length_set(sequence) == lengths
        for k in range(len(terms) // 2 + 2):
            assert local.has_factorization_of_length(k) == (k in lengths)
            witness = local.factorization_witness(k)
            assert (witness is not None) == (k in lengths)
            if witness is not None:
                assert sum(witness, space()) == sequence
        found = {tuple(sorted(atom_ids[a] for a in f)) for f in local.factorizations()}
        assert found == set(expected)
        if len(terms) <= 3:
            graph = local.digraph()
            assert nx.is_directed_acyclic_graph(graph)
            assert all(t + data["atom"] == s for s, t, data in graph.edges(data=True))


@pytest.mark.parametrize("moduli", GROUPS, ids=str)
def test_orbits_and_canonical_forms(moduli):
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(*moduli))
    group = space.automorphism_group()
    inputs = list(multisets(elements_of(moduli), 4))
    for terms in inputs[:: max(1, len(inputs) // 9)]:
        sequence = space(terms)
        images = {group.apply(g, sequence) for g in group}
        assert set(sequence.orbit()) == images
        form = group.canonical_form(sequence)
        assert form.image == min(images) == group.apply(form.automorphism, sequence)


def cyclic_subgroup(group, generator):
    powers, current = [group.identity], generator
    while current != group.identity:
        powers.append(current)
        current = group.compose(generator, current)
    return group.subgroup(powers)


@pytest.mark.parametrize("moduli", [(1,), (2,), (3,), (4,), (2, 2)], ids=str)
@pytest.mark.parametrize("lengths", [
    {1}, {2}, {3}, {4}, {1, 2}, {1, 3}, {2, 3}, {2, 4}, {1, 2, 3, 4},
], ids=str)
def test_zero_sum_constants(moduli, lengths):
    elements = elements_of(moduli)
    zero = (0,) * len(moduli)
    multiples = [k for k in lengths if k % lcm(*moduli) == 0]

    def avoids(terms):
        return not any(
            len(sub) in lengths and vector_sum(sub, moduli) == zero
            for sub in sub_multisets(terms)
        )

    value, extremal = inf, ()
    if multiples:
        # Pigeonhole: this many terms contain min(multiples) equal ones.
        previous = [()]
        for k in range(1, len(elements) * (min(multiples) - 1) + 2):
            free = [t for t in combinations_with_replacement(elements, k) if avoids(t)]
            if not free:
                value, extremal = k, previous
                break
            previous = free

    # An insufficient explicit bound must not influence the constants.
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(*moduli), davenport_bound=1)
    full = space.automorphism_group()
    subgroup = cyclic_subgroup(full, full.generators[0]) if full.generators else full
    all_extremal = {space(t) for t in extremal}
    for symmetry, group in ((False, None), (True, full), (subgroup, subgroup)):
        result = zero_sum_constant(space, lengths, symmetry=symmetry)
        expected = all_extremal if group is None else {group.canonical(s) for s in all_extremal}
        assert result.value == value
        assert set(result.extremal) == expected
        assert space.davenport_bound == 1
