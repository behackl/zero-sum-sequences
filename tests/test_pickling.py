"""Pickling keeps sequences usable across processes."""

import copy
import multiprocessing
import pickle
from operator import xor

from zero_sum_sequences import (
    AdditiveSequenceSpace,
    FactorizationRelation,
    FiniteAdditiveGroup,
    davenport_constant,
)


def _length_set(sequence):
    return sequence, sequence.length_set(), sequence.parent().davenport_bound


def _catalogue(space):
    return space.enumerate_atom_catalogue()


def _space():
    return AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(2, 4))


def test_round_trips_in_one_process_preserve_identity():
    space = _space()
    sequence = space([(0, 1), (0, 3)])
    assert pickle.loads(pickle.dumps(space)) is space
    assert pickle.loads(pickle.dumps(space.base_parent)) is space.base_parent
    assert pickle.loads(pickle.dumps(sequence)) == sequence
    assert copy.deepcopy(sequence) == sequence


def test_unpickled_copies_are_rebuilt_from_the_recipe():
    space = _space()
    davenport_constant(space)
    data = pickle.dumps(space([(0, 1), (0, 3)]))
    del space  # no live original: unpickling builds a new space
    sequence = pickle.loads(data)
    assert sequence.is_atom() and sequence.parent().davenport_bound == 5
    assert sequence.parent().base_parent._group_arithmetic().exact_davenport_constant == 5


def test_custom_groups_pickle_when_their_callables_do():
    space = AdditiveSequenceSpace(FiniteAdditiveGroup(range(4), zero=0, add=xor))
    restored = pickle.loads(pickle.dumps(space([1, 1])))
    assert restored.is_atom() and restored.parent().base_parent.exponent() == 2


def test_catalogues_groups_and_relations_round_trip():
    space = _space()
    catalogue = space.enumerate_atom_catalogue()
    assert pickle.loads(pickle.dumps(catalogue)).digest() == catalogue.digest()
    group = space.automorphism_group()
    assert pickle.loads(pickle.dumps(group)).elements == group.elements
    product = catalogue[-1] + catalogue[-2]
    factorizations = list(product.factorizations())
    relation = FactorizationRelation(factorizations[0], factorizations[-1])
    assert pickle.loads(pickle.dumps(relation)) == relation


def test_worker_results_belong_to_the_callers_space():
    space = _space()
    inputs = [space([(0, 1)] * 4), space([(0, 1), (0, 3)] * 2)]
    with multiprocessing.get_context("spawn").Pool(2) as pool:
        results = pool.map(_length_set, inputs)
        (catalogue,) = pool.map(_catalogue, [space])
    for (sequence, lengths, bound), original in zip(results, inputs):
        assert sequence.parent() is space and sequence == original
        assert lengths == original.length_set() and bound == 5
    assert catalogue.space is space
    assert catalogue.atoms == space.enumerate_atom_catalogue().atoms
