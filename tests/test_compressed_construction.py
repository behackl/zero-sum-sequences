"""Multiplicity construction stays compressed without weakening canonicalization."""

import pytest

from zero_sum_sequences import AdditiveSequence, AdditiveSequenceSpace, FiniteAdditiveGroup


def test_large_multiplicity_coerces_the_key_once():
    class Parent:
        calls = 0

        def __call__(self, term):
            self.calls += 1
            return int(term)

        def zero(self):
            return 0

    parent = Parent()
    space = AdditiveSequenceSpace(parent)
    sequence = space.from_multiplicities({1: 100_000})
    assert parent.calls == 1
    assert sequence._items == ((1, 100_000),)
    assert len(sequence) == 100_000
    assert space._arithmetic_context is None


def test_compressed_construction_does_not_expand_via_the_general_constructor(monkeypatch):
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(3))

    def forbidden(*args):
        raise AssertionError("multiplicity construction must not expand terms")

    monkeypatch.setattr(AdditiveSequence, "__init__", forbidden)
    sequence = space.from_multiplicities({(1,): 100_000})
    assert sequence._items == (((1,), 100_000),)
    assert len(sequence) == 100_000


def test_aliases_are_merged_and_zero_count_keys_are_not_coerced():
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(3))
    result = space.from_multiplicities({(1,): 2, (4,): 3, (2,): 1, object(): 0})
    expected = space([(1,)] * 5 + [(2,)])
    assert result == expected
    assert hash(result) == hash(expected)
    assert result.multiplicities == {(1,): 5, (2,): 1}
    assert result.parent() is space
    assert space.from_multiplicities({object(): 0}) == space()


@pytest.mark.parametrize("bad", [True, -1, 1.5])
def test_compressed_construction_rejects_invalid_counts(bad):
    space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(3))
    with pytest.raises(ValueError, match="multiplicities"):
        space.from_multiplicities({(1,): bad})


def test_nonorderable_terms_keep_the_construction_error():
    class Parent:
        def __call__(self, term):
            return term

        def zero(self):
            return 0

    with pytest.raises(TypeError, match="mutually orderable"):
        AdditiveSequenceSpace(Parent()).from_multiplicities({1: 1, "a": 1})
