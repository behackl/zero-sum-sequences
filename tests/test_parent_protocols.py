"""Foreign-parent arithmetic and capability failure boundaries."""

from functools import reduce

import pytest

from zero_sum_sequences import (
    AdditiveSequenceSpace,
    AtomCatalogue,
    AutomorphismGroup,
    AutomorphismGroupUnavailable,
    FiniteAdditiveGroup,
    davenport_constant,
    zero_sum_constant,
)


class Label(int):
    def is_zero(self):
        return int(self) == 0


class ForeignC3:
    def __call__(self, term):
        return int(term) % 3

    def zero(self):
        return 0

    def is_finite(self):
        return True

    def __iter__(self):
        return iter(range(3))

    def add(self, left, right):
        assert 0 <= left < 3 and 0 <= right < 3
        return left + right  # The parent canonicalizes the result.


@pytest.mark.parametrize("parent", [
    lambda: FiniteAdditiveGroup(
        map(Label, [1, 2, 3]), zero=Label(1),
        coerce=lambda x: Label((int(x) - 1) % 3 + 1),
        add=lambda a, b: Label((int(a) + int(b) - 2) % 3 + 1),
    ),
    ForeignC3,
])
def test_parent_arithmetic_agrees_across_predicates_catalogues_and_constants(parent):
    parent = parent()
    space = AdditiveSequenceSpace(parent)
    nonzero = tuple(t for t in parent if t != parent.zero())
    a, b = nonzero
    pair, pure_a, pure_b = space(a, b), space(a) ** 3, space(b) ** 3
    for sequence in (pair, pure_a, pure_b):
        assert sequence.sum() == parent.zero()
        assert sequence.is_zero_sum()
        assert sequence.is_atom()
    assert not space(a).is_zero_sum()
    assert not (pair * pair).is_atom()
    assert space().is_zero_sum()
    catalogue = space.enumerate_atom_catalogue()
    assert set(catalogue) == {pair, pure_a, pure_b}
    sequence = pure_a * pure_b
    assert sequence.length_set() == sequence.length_set(atom_catalogue=catalogue) == {2, 3}
    assert davenport_constant(space, symmetry=False) == 3
    assert zero_sum_constant(space, {1, 2, 3}, symmetry=False).value == 3


def test_identity_is_the_coerced_zero_label():
    class ShiftedZero(ForeignC3):
        def zero(self):
            return 3  # A non-canonical label of the identity 0.

    space = AdditiveSequenceSpace(ShiftedZero())
    assert space().sum() == 0
    assert space(1, 2).is_zero_sum() and space(1, 2).is_atom()
    assert len(space.enumerate_atom_catalogue()) == 3
    assert davenport_constant(space, symmetry=False) == 3


def test_unknown_finiteness_is_handled_alike_everywhere():
    class UnknownFiniteness(ForeignC3):
        def is_finite(self):
            raise NotImplementedError

    space = AdditiveSequenceSpace(UnknownFiniteness())
    with pytest.raises(ValueError, match="supply davenport_bound explicitly"):
        _ = space.davenport_bound
    with pytest.raises(ValueError, match="requires a finite parent"):
        space.enumerate_atom_catalogue()
    assert space(1, 2).length_set() == {1}


def test_native_group_totals_do_not_coerce_an_already_canonical_result():
    coercions = []

    def coerce(term):
        coercions.append(term)
        return int(term) % 3

    parent = FiniteAdditiveGroup(range(3), zero=0, add=lambda a, b: (a + b) % 3, coerce=coerce)
    space = AdditiveSequenceSpace(parent, davenport_bound=3)
    sequence = space(1, 2)
    zero = space().sum()
    coercions.clear()
    expected = reduce(parent.add, sequence, zero)
    native_coercions = len(coercions)
    coercions.clear()
    assert sequence.sum() == expected
    assert len(coercions) == native_coercions


def test_totals_do_not_require_a_finite_parent():
    class Integers:
        def __call__(self, term):
            return int(term)

        def zero(self):
            return 0

    space = AdditiveSequenceSpace(Integers())
    assert space(-2, 1, 1).is_atom()


class MetadataC3:
    def __call__(self, term):
        return int(term) % 3

    def zero(self):
        return 0

    def add(self, left, right):
        return (left + right) % 3

    def is_finite(self):
        return True

    def __iter__(self):
        return iter(range(3))


def parent_with_failure(provider, error):
    def broken(*args):
        raise error

    return type("BrokenParent", (MetadataC3,), {provider: broken})()


@pytest.mark.parametrize(
    "provider, error_type",
    [
        ("cardinality", ValueError),
        ("order", ValueError),
        ("invariants", ValueError),
        ("exponent", ValueError),
        ("is_finite", ValueError),
        ("__len__", TypeError),
        ("__iter__", TypeError),
        ("cardinality", AttributeError),
        ("is_finite", RuntimeError),  # NotImplementedError is a RuntimeError subclass.
    ],
)
def test_bound_provider_errors_are_not_unavailability(provider, error_type):
    error = error_type("broken metadata provider")
    space = AdditiveSequenceSpace(parent_with_failure(provider, error))
    with pytest.raises(error_type) as raised:
        _ = space.davenport_bound
    assert raised.value is error


@pytest.mark.parametrize("query", ["factorization", "enumeration", "export"])
def test_bound_consumers_do_not_hide_broken_metadata(query, tmp_path):
    error = ValueError("broken metadata provider")
    space = AdditiveSequenceSpace(parent_with_failure("cardinality", error))
    sequence = space(1, 2)
    catalogue = AtomCatalogue(space, [sequence])
    queries = {
        "factorization": sequence.factorization_solver,
        "enumeration": space.enumerate_atom_catalogue,
        "export": lambda: catalogue.to_jsonl(tmp_path / "catalogue.jsonl"),
    }
    with pytest.raises(ValueError) as raised:
        queries[query]()
    assert raised.value is error


@pytest.mark.parametrize("metadata, message", [
    ({"cardinality": lambda self: -1}, "group order"),
    ({"invariants": lambda self: (2,)}, "do not match"),
    ({"exponent": lambda self: 2}, "divide"),
])
def test_invalid_metadata_is_not_silently_ignored(metadata, message):
    space = AdditiveSequenceSpace(type("InvalidParent", (MetadataC3,), metadata)())
    with pytest.raises(ValueError, match=message):
        space(1, 2).factorization_solver()


def test_unsupported_metadata_uses_remaining_capabilities():
    def unsupported(*args):
        raise NotImplementedError

    cls = type("PartialParent", (MetadataC3,), {
        "cardinality": unsupported, "order": unsupported,
        "__len__": unsupported, "invariants": unsupported, "exponent": unsupported,
    })
    space = AdditiveSequenceSpace(cls())
    assert space.davenport_bound == 3
    assert space(1, 2).length_set() == {1}


@pytest.mark.parametrize("iterator", [None, lambda self: (_ for _ in ()).throw(NotImplementedError())])
def test_explicitly_unsupported_iteration_allows_local_fallback(iterator):
    cls = type("NoIteration", (MetadataC3,), {"__iter__": iterator})
    space = AdditiveSequenceSpace(cls())
    assert space._optional_davenport_bound() is None
    assert space(1, 2).length_set() == {1}


@pytest.mark.parametrize("provider", ["__iter__", "__len__"])
def test_invalid_iteration_provider_does_not_allow_local_fallback(provider):
    space = AdditiveSequenceSpace(type("InvalidIteration", (MetadataC3,), {provider: 5})())
    with pytest.raises(TypeError):
        space(1, 2).factorization_solver()


def test_sequence_indexing_protocol_still_supplies_iteration():
    class Indexing:
        __call__ = MetadataC3.__call__
        zero = MetadataC3.zero
        add = MetadataC3.add
        is_finite = MetadataC3.is_finite

        def __getitem__(self, index):
            return range(3)[index]

    space = AdditiveSequenceSpace(Indexing())
    assert space.davenport_bound == 3
    assert space(1, 2).length_set() == {1}


def test_missing_metadata_still_allows_boundless_local_solving():
    class Limited:
        def __call__(self, value):
            return int(value) % 3

        def zero(self):
            return 0

        def add(self, a, b):
            return (a + b) % 3

        def is_finite(self):
            return True

    space = AdditiveSequenceSpace(Limited())
    assert space._optional_davenport_bound() is None
    assert space(1, 2).length_set() == {1}
    with pytest.raises(ValueError, match="supply davenport_bound explicitly"):
        _ = space.davenport_bound


@pytest.mark.parametrize(
    "provider, error_type",
    [
        ("automorphism_group", AttributeError),
        ("automorphism_action", AttributeError),
        ("automorphism_generators", AttributeError),
        ("automorphism_group", TypeError),
        ("automorphism_action", ValueError),
        ("automorphism_generators", RuntimeError),
    ],
)
def test_explicit_symmetry_provider_failures_propagate(provider, error_type):
    error = error_type("broken symmetry provider")

    space = AdditiveSequenceSpace(parent_with_failure(provider, error))
    with pytest.raises(error_type) as raised:
        zero_sum_constant(space, {3})
    assert raised.value is error
    assert zero_sum_constant(space, {3}, symmetry=False).value == 5


def test_image_table_construction_preserves_iterator_provider_errors():
    error = TypeError("broken iterator")
    parent = parent_with_failure("__iter__", error)
    with pytest.raises(TypeError) as raised:
        AutomorphismGroup.closure(parent, [lambda x: x])
    assert raised.value is error


def test_configured_package_group_provider_errors_propagate():
    def broken():
        raise AttributeError("configured callback failed")

    parent = FiniteAdditiveGroup(range(3), zero=0, add=lambda a, b: (a + b) % 3,
                                 automorphism_group=broken)
    with pytest.raises(AttributeError, match="configured callback"):
        zero_sum_constant(AdditiveSequenceSpace(parent), {3})


def test_unsupported_symmetry_still_falls_back():
    def unsupported(*args):
        raise NotImplementedError

    cls = type("UnsupportedSymmetry", (MetadataC3,), {
        name: unsupported for name in (
            "automorphism_group", "automorphism_action", "automorphism_generators"
        )
    })
    assert zero_sum_constant(AdditiveSequenceSpace(cls()), {3}).value == 5


def _negation(term):
    return (-term) % 3


class MissingIterationC3:
    __call__ = MetadataC3.__call__
    zero = MetadataC3.zero
    add = MetadataC3.add
    is_finite = MetadataC3.is_finite

    def cardinality(self):
        return 3

    def automorphism_generators(self):
        return (_negation,)


@pytest.mark.parametrize(
    "parent_type",
    [
        MissingIterationC3,
        type("DisabledIteration", (MissingIterationC3,), {"__iter__": None}),
        type(
            "UnsupportedIteration",
            (MissingIterationC3,),
            {"__iter__": lambda self: (_ for _ in ()).throw(NotImplementedError())},
        ),
    ],
    ids=["missing", "disabled", "not-implemented"],
)
def test_symmetry_falls_back_when_image_tables_need_iteration(parent_type):
    space = AdditiveSequenceSpace(parent_type())
    a, b, c = space(1, 1, 1), space(2, 2, 2), space(1, 2)
    catalogue = AtomCatalogue(space, [a, b, c])

    assert set(a.orbit()) == {a, b}
    with pytest.raises(AutomorphismGroupUnavailable):
        space.automorphism_group()
    cache = space.factorization_cache(catalogue, symmetry=True)
    assert cache.group is None
    assert cache.length_set(a * b) == {2, 3}


def test_broken_iteration_is_not_mistaken_for_unavailable_symmetry():
    error = TypeError("broken iterator")

    def broken(self):
        raise error

    space = AdditiveSequenceSpace(
        type("BrokenIteration", (MissingIterationC3,), {"__iter__": broken})()
    )
    catalogue = AtomCatalogue(space, [space(1, 1, 1), space(2, 2, 2), space(1, 2)])
    with pytest.raises(TypeError) as raised:
        space.factorization_cache(catalogue, symmetry=True)
    assert raised.value is error
