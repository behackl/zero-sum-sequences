"""Shared short-to-very-long factorization benchmark corpus."""

from __future__ import annotations

from dataclasses import dataclass

from sage.all import GF, Zmod

from zero_sum_sequences import AdditiveSequence, AdditiveSequenceSpace


@dataclass(frozen=True, slots=True)
class FactorizationBenchmarkCase:
    name: str
    tier: str
    sequence: AdditiveSequence
    expected_lengths: frozenset[int]
    enumerate_factorizations: bool = False
    expected_factorizations: int | None = None


def factorization_benchmark_cases() -> tuple[FactorizationBenchmarkCase, ...]:
    """Return deterministic benchmark inputs with mathematical expectations."""

    c3 = Zmod(3)
    c3_sequences = AdditiveSequenceSpace(c3, davenport_bound=3)
    positive_atom = c3_sequences([c3(1)] * 3)
    negative_atom = c3_sequences([c3(2)] * 3)
    balanced_block = positive_atom + negative_atom

    group = GF(3) ** 3
    c3_cubed_sequences = AdditiveSequenceSpace(group, davenport_bound=7)
    e1, e2, e3 = group.basis()
    rank_three_atom = c3_cubed_sequences(
        [e1, e1, e2, e2, e3, e3, e1 + e2 + e3]
    )
    negative_rank_three_atom = c3_cubed_sequences(
        -term for term in rank_three_atom
    )
    rank_three_inverse_pair = rank_three_atom + negative_rank_three_atom

    return (
        FactorizationBenchmarkCase(
            name="c3-atom",
            tier="short",
            sequence=positive_atom,
            expected_lengths=frozenset({1}),
            enumerate_factorizations=True,
            expected_factorizations=1,
        ),
        FactorizationBenchmarkCase(
            name="c3-balanced-block",
            tier="short",
            sequence=balanced_block,
            expected_lengths=frozenset({2, 3}),
            enumerate_factorizations=True,
            expected_factorizations=2,
        ),
        FactorizationBenchmarkCase(
            name="rank-three-inverse-pair",
            tier="short",
            sequence=rank_three_inverse_pair,
            expected_lengths=frozenset({2, 3, 4, 5, 7}),
            enumerate_factorizations=False,
        ),
        FactorizationBenchmarkCase(
            name="c3-balanced-power-7",
            tier="long",
            sequence=7 * balanced_block,
            expected_lengths=frozenset(range(14, 22)),
            enumerate_factorizations=True,
            expected_factorizations=8,
        ),
        FactorizationBenchmarkCase(
            name="c3-balanced-power-40",
            tier="long",
            sequence=40 * balanced_block,
            expected_lengths=frozenset(range(80, 121)),
            enumerate_factorizations=True,
            expected_factorizations=41,
        ),
        FactorizationBenchmarkCase(
            name="rank-three-inverse-power-3",
            tier="long",
            sequence=3 * rank_three_inverse_pair,
            expected_lengths=frozenset(range(6, 22)),
        ),
        FactorizationBenchmarkCase(
            name="c3-balanced-power-250",
            tier="very-long",
            sequence=250 * balanced_block,
            expected_lengths=frozenset(range(500, 751)),
            enumerate_factorizations=True,
            expected_factorizations=251,
        ),
        FactorizationBenchmarkCase(
            name="c3-pure-power-5000",
            tier="very-long",
            sequence=5000 * positive_atom,
            expected_lengths=frozenset({5000}),
            enumerate_factorizations=True,
            expected_factorizations=1,
        ),
    )
