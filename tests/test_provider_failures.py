"""Unavailable capabilities permit fallback; broken providers do not."""

import pytest

from zero_sum_sequences import (
    AdditiveSequenceSpace,
    AtomCatalogue,
    AutomorphismGroup,
    FiniteAdditiveGroup,
    zero_sum_constant,
)


class Parent:
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


@pytest.mark.parametrize("provider", ["cardinality", "order", "invariants", "exponent", "is_finite", "__len__", "__iter__"])
@pytest.mark.parametrize("error_type", [TypeError, ValueError, AttributeError, RuntimeError])
def test_bound_provider_errors_are_not_unavailability(provider, error_type, tmp_path):
    error = error_type("broken metadata provider")

    def broken(*args):
        raise error

    cls = type("BrokenParent", (Parent,), {provider: broken})
    space = AdditiveSequenceSpace(cls())
    sequence = space([1, 2])
    catalogue = AtomCatalogue(space, [sequence])
    for query in (
        lambda: space.davenport_bound,
        space._optional_davenport_bound,
        sequence.factorization_solver,
        space.enumerate_atom_catalogue,
        lambda: catalogue.to_jsonl(tmp_path / "catalogue.jsonl"),
    ):
        with pytest.raises(error_type) as raised:
            query()
        assert raised.value is error


@pytest.mark.parametrize("metadata, message", [
    ({"cardinality": lambda self: -1}, "group order"),
    ({"invariants": lambda self: (2,)}, "do not match"),
    ({"exponent": lambda self: 2}, "divide"),
])
def test_invalid_metadata_is_not_silently_ignored(metadata, message):
    space = AdditiveSequenceSpace(type("InvalidParent", (Parent,), metadata)())
    with pytest.raises(ValueError, match=message):
        space([1, 2]).factorization_solver()


def test_unsupported_metadata_uses_remaining_capabilities():
    def unsupported(*args):
        raise NotImplementedError

    cls = type("PartialParent", (Parent,), {
        "cardinality": unsupported, "order": unsupported,
        "__len__": unsupported, "invariants": unsupported, "exponent": unsupported,
    })
    space = AdditiveSequenceSpace(cls())
    assert space.davenport_bound == 3
    assert space([1, 2]).length_set() == {1}


@pytest.mark.parametrize("iterator", [None, lambda self: (_ for _ in ()).throw(NotImplementedError())])
def test_explicitly_unsupported_iteration_allows_local_fallback(iterator):
    cls = type("NoIteration", (Parent,), {"__iter__": iterator})
    space = AdditiveSequenceSpace(cls())
    assert space._optional_davenport_bound() is None
    assert space([1, 2]).length_set() == {1}


@pytest.mark.parametrize("provider", ["__iter__", "__len__"])
def test_invalid_iteration_provider_does_not_allow_local_fallback(provider):
    space = AdditiveSequenceSpace(type("InvalidIteration", (Parent,), {provider: 5})())
    with pytest.raises(TypeError):
        space([1, 2]).factorization_solver()


def test_sequence_indexing_protocol_still_supplies_iteration():
    class Indexing:
        __call__ = Parent.__call__
        zero = Parent.zero
        add = Parent.add
        is_finite = Parent.is_finite

        def __getitem__(self, index):
            return range(3)[index]

    space = AdditiveSequenceSpace(Indexing())
    assert space.davenport_bound == 3
    assert space([1, 2]).length_set() == {1}


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
    assert space([1, 2]).length_set() == {1}
    with pytest.raises(ValueError, match="supply davenport_bound explicitly"):
        _ = space.davenport_bound


@pytest.mark.parametrize("provider", ["automorphism_group", "automorphism_action", "automorphism_generators"])
@pytest.mark.parametrize("error_type", [AttributeError, TypeError, ValueError, RuntimeError])
def test_explicit_symmetry_provider_failures_propagate(provider, error_type):
    error = error_type("broken symmetry provider")

    def broken(*args):
        raise error

    space = AdditiveSequenceSpace(type("BrokenSymmetry", (Parent,), {provider: broken})())
    with pytest.raises(error_type) as raised:
        zero_sum_constant(space, {3})
    assert raised.value is error
    assert zero_sum_constant(space, {3}, symmetry=False).value == 5


@pytest.mark.parametrize("error_type", [AttributeError, TypeError, ValueError, RuntimeError])
def test_image_table_construction_preserves_iterator_provider_errors(error_type):
    error = error_type("broken iterator")

    def broken(self):
        raise error

    parent = type("BrokenIterator", (Parent,), {"__iter__": broken})()
    with pytest.raises(error_type) as raised:
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

    cls = type("UnsupportedSymmetry", (Parent,), {
        name: unsupported for name in (
            "automorphism_group", "automorphism_action", "automorphism_generators"
        )
    })
    assert zero_sum_constant(AdditiveSequenceSpace(cls()), {3}).value == 5
