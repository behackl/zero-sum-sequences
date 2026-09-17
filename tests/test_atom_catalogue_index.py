"""Indexing, automorphism orbits and serialization of atom catalogues."""

from __future__ import annotations

import json

import pytest

from zero_sum_sequences import (
    AdditiveSequenceSpace,
    AtomCatalogue,
    AtomOrbit,
    FiniteAdditiveGroup,
)


def space_for(*moduli, bound):
    return AdditiveSequenceSpace(
        FiniteAdditiveGroup.cyclic_product(*moduli), davenport_bound=bound
    )


@pytest.fixture(params=[((2, 4), 5), ((6,), 6), ((3, 3), 5), ((2, 2, 2), 4)])
def catalogue(request):
    moduli, bound = request.param
    return space_for(*moduli, bound=bound).enumerate_atom_catalogue()


# order and indexing


def test_catalogue_order_is_the_sequence_order(catalogue):
    atoms = catalogue.atoms
    assert list(atoms) == sorted(atoms)
    for position, atom in enumerate(atoms):
        assert catalogue.index(atom) == position
        assert catalogue[position] is atom
        assert atom in catalogue
    not_an_atom = catalogue.space(tuple(atoms[0]) + tuple(atoms[1]))
    assert not_an_atom not in catalogue
    with pytest.raises(ValueError):
        catalogue.index(not_an_atom)


def test_filter_and_restrict(catalogue):
    short = catalogue.restrict(max_length=2)
    assert all(len(atom) <= 2 for atom in short)
    assert len(short) == sum(len(atom) <= 2 for atom in catalogue)
    first = catalogue[-1]
    on_support = catalogue.restrict(support=first.support)
    assert first in on_support
    assert all(set(atom.support) <= set(first.support) for atom in on_support)
    assert catalogue.filter(lambda atom: False).atoms == ()
    assert list(catalogue.restrict(max_length=2).divisors(first)) == [
        atom for atom in catalogue.divisors(first) if len(atom) <= 2
    ]


# automorphisms


def test_permutation_is_a_permutation_compatible_with_the_action(catalogue):
    group = catalogue.space.automorphism_group()
    for element in group.elements[:5]:
        permutation = catalogue.permutation(element)
        assert sorted(permutation) == list(range(len(catalogue)))
        for position, atom in enumerate(catalogue):
            assert catalogue[permutation[position]] == group.apply(element, atom)
    assert catalogue.permutation(group.identity) == tuple(range(len(catalogue)))
    # composition of permutations matches composition of automorphisms
    a, b = group.elements[1], group.elements[-1]
    pa, pb, pab = (catalogue.permutation(x) for x in (a, b, group.compose(a, b)))
    assert pab == tuple(pa[pb[i]] for i in range(len(catalogue)))


def test_permutation_requires_closure():
    catalogue = space_for(2, 4, bound=5).enumerate_atom_catalogue()
    group = catalogue.space.automorphism_group()
    partial = catalogue.filter(lambda atom: atom == catalogue[-1])
    with pytest.raises(ValueError):
        partial.permutation(group.elements[1])


@pytest.mark.parametrize(
    "moduli, bound, orbits",
    [((2, 4), 5, 11), ((3, 3), 5, 5), ((2, 2, 2), 4, 3), ((6,), 6, 11), ((2, 2, 2, 2), 5, 4)],
)
def test_orbits(moduli, bound, orbits):
    catalogue = space_for(*moduli, bound=bound).enumerate_atom_catalogue()
    group = catalogue.space.automorphism_group()
    result = catalogue.orbits()
    assert len(result) == orbits
    assert all(isinstance(orbit, AtomOrbit) for orbit in result)
    covered = sorted(index for orbit in result for index in orbit.indices)
    assert covered == list(range(len(catalogue)))
    for orbit in result:
        assert len(orbit) * orbit.stabilizer_order == len(group)
        assert orbit.representative == min(catalogue[i] for i in orbit.indices)
        assert all(
            group.canonical(catalogue[i]) == orbit.representative for i in orbit.indices
        )
    assert [orbit.representative for orbit in result] == sorted(
        orbit.representative for orbit in result
    )


def test_canonical_delegates_to_the_group():
    catalogue = space_for(2, 4, bound=5).enumerate_atom_catalogue()
    atom = catalogue[-1]
    form = catalogue.canonical(atom)
    assert form.image in catalogue
    assert catalogue.space.automorphism_group().apply(form.automorphism, atom) == form.image


# serialization


def test_jsonl_round_trip_with_annotations(catalogue, tmp_path):
    path = tmp_path / "catalogue.jsonl"
    group = catalogue.space.automorphism_group()
    written = catalogue.to_jsonl(
        path, annotate=lambda atom: {"orbit_size": len(group.orbit(atom))}
    )
    assert written == len(catalogue)
    lines = path.read_text().splitlines()
    metadata = json.loads(lines[0])
    assert metadata["count"] == len(catalogue) and metadata["digest"] == catalogue.digest()
    records = [json.loads(line) for line in lines[1:]]
    assert [r["index"] for r in records] == list(range(len(catalogue)))
    assert all(r["length"] == len(r["sequence"]) for r in records)

    loaded = AtomCatalogue.from_jsonl(path, catalogue.space)
    assert loaded.atoms == catalogue.atoms
    assert loaded.digest() == catalogue.digest()
    assert all(
        loaded.annotations[atom]["orbit_size"] == len(group.orbit(atom)) for atom in loaded
    )
    fast = AtomCatalogue.from_jsonl(path, catalogue.space, verify=False)
    assert fast.atoms == loaded.atoms and fast.annotations == loaded.annotations
    assert list(fast.divisors(catalogue[-1])) == list(catalogue.divisors(catalogue[-1]))
    # annotations survive a round trip through records() again
    assert next(iter(loaded.records()))["orbit_size"] == next(iter(catalogue.records(
        annotate=lambda atom: {"orbit_size": len(group.orbit(atom))}
    )))["orbit_size"]


def test_digest_ignores_annotations_and_depends_on_atoms():
    space = space_for(2, 4, bound=5)
    catalogue = space.enumerate_atom_catalogue()
    annotated = AtomCatalogue(space, catalogue.atoms, annotations=[(catalogue[0], {"a": 1})])
    assert annotated.digest() == catalogue.digest()
    assert catalogue.restrict(max_length=3).digest() != catalogue.digest()
    assert len(catalogue.digest()) == 64


def test_from_jsonl_rejects_tampering(tmp_path):
    space = space_for(2, 4, bound=5)
    catalogue = space.enumerate_atom_catalogue()
    path = tmp_path / "catalogue.jsonl"
    catalogue.to_jsonl(path)
    lines = path.read_text().splitlines()

    def write(modified):
        path.write_text("\n".join(modified) + "\n")

    # swapped records: order violation
    swapped = [lines[0], lines[2], lines[1], *lines[3:]]
    write(swapped)
    with pytest.raises(ValueError):
        AtomCatalogue.from_jsonl(path, space)
    # a non-atom sequence with verify=True
    record = json.loads(lines[1])
    record["sequence"] = record["sequence"] + record["sequence"]
    record["length"] = len(record["sequence"])
    write([lines[0], json.dumps(record), *lines[2:]])
    with pytest.raises(ValueError):
        AtomCatalogue.from_jsonl(path, space)
    # dropped record: count and digest mismatch
    write(lines[:-1])
    with pytest.raises(ValueError):
        AtomCatalogue.from_jsonl(path, space)
    # wrong schema
    write(['{"schema":"other"}', *lines[1:]])
    with pytest.raises(ValueError):
        AtomCatalogue.from_jsonl(path, space)
    # reserved annotation key
    with pytest.raises(ValueError):
        list(catalogue.records(annotate=lambda atom: {"index": 0}))


def test_annotations_must_belong_to_the_catalogue():
    space = space_for(6, bound=6)
    catalogue = space.enumerate_atom_catalogue()
    stranger = space([(1,), (5,), (2,), (4,)])
    with pytest.raises(ValueError):
        AtomCatalogue(space, catalogue.atoms, annotations=[(stranger, {})])


def test_term_codec_on_parents():
    parent = FiniteAdditiveGroup.cyclic_product(2, 4)
    assert parent.encode_term((1, 3)) == [1, 3]
    assert parent.decode_term([1, 7]) == (1, 3)
    custom = FiniteAdditiveGroup(
        range(5),
        zero=0,
        add=lambda a, b: (a + b) % 5,
        encode_term=lambda t: f"t{t}",
        decode_term=lambda d: int(d[1:]),
    )
    assert custom.encode_term(3) == "t3" and custom.decode_term("t3") == 3
    space = AdditiveSequenceSpace(custom, davenport_bound=5)
    x = space([1, 1, 3])
    assert x.encode() == ["t1", "t1", "t3"] and space.decode(x.encode()) == x
    with pytest.raises(TypeError):
        FiniteAdditiveGroup(range(5), zero=0, add=lambda a, b: (a + b) % 5, encode_term=str)
