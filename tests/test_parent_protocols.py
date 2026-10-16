"""Cross-layer arithmetic contracts for ordinary and foreign parent labels."""

import pytest

from zero_sum_sequences import (
    AdditiveSequenceSpace,
    AtomCatalogue,
    AutomorphismGroup,
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
    pair, pure_a, pure_b = space([a, b]), space([a] * 3), space([b] * 3)
    for sequence in (pair, pure_a, pure_b):
        assert sequence.total() == parent.zero()
        assert sequence.is_zero_sum()
        assert sequence.is_atom()
    assert not space([a]).is_zero_sum()
    assert not (pair + pair).is_atom()
    assert space().is_zero_sum()
    assert space._arithmetic_context is None
    catalogue = space.enumerate_atom_catalogue()
    assert set(catalogue) == {pair, pure_a, pure_b}
    sequence = pure_a + pure_b
    assert sequence.length_set() == sequence.length_set(atom_catalogue=catalogue) == {2, 3}
    assert davenport_constant(space, symmetry=False) == 3
    assert zero_sum_constant(space, {1, 2, 3}, symmetry=False).value == 3


def test_identity_is_the_coerced_zero_label():
    class ShiftedZero(ForeignC3):
        def zero(self):
            return 3  # A non-canonical label of the identity 0.

    space = AdditiveSequenceSpace(ShiftedZero())
    assert space().total() == 0
    assert space([1, 2]).is_zero_sum() and space([1, 2]).is_atom()
    assert len(space.enumerate_atom_catalogue()) == 3
    assert davenport_constant(space, symmetry=False) == 3


def test_unknown_finiteness_is_handled_alike_everywhere():
    class UnknownFiniteness(ForeignC3):
        def is_finite(self):
            raise NotImplementedError

    space = AdditiveSequenceSpace(UnknownFiniteness())
    assert space._optional_davenport_bound() is None
    with pytest.raises(ValueError, match="requires a finite parent"):
        space.enumerate_atom_catalogue()
    assert space([1, 2]).length_set() == {1}


def test_native_group_totals_do_not_coerce_an_already_canonical_result():
    coercions = []

    def coerce(term):
        coercions.append(term)
        return int(term) % 3

    parent = FiniteAdditiveGroup(range(3), zero=0, add=lambda a, b: (a + b) % 3, coerce=coerce)
    space = AdditiveSequenceSpace(parent, davenport_bound=3)
    sequence = space([1, 2])
    space().total()  # Resolves the canonical identity once.
    coercions.clear()
    assert sequence.total() == 0
    assert len(coercions) == 3 * len(sequence)  # FiniteAdditiveGroup.add's own checks only.
    assert parent._arithmetic_context is space._arithmetic_context is None


def test_totals_do_not_require_a_finite_parent():
    class Integers:
        def __call__(self, term):
            return int(term)

        def zero(self):
            return 0

    space = AdditiveSequenceSpace(Integers())
    assert space([-2, 1, 1]).is_atom()
    assert space._arithmetic_context is None


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
