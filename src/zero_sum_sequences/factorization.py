"""Memoized factorizations of additive sequences into atoms."""

from __future__ import annotations

from collections import deque
from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass
from typing import Generic, TypeAlias

import networkx as nx

from .additive_sequence import AdditiveSequence, Element, _non_negative_integer
from .atom_catalogue import AtomCatalogue

Factorization: TypeAlias = tuple[AdditiveSequence[Element], ...]
MultiplicityState: TypeAlias = tuple[int, ...]


def _atom_divisors(
    sequence: AdditiveSequence[Element],
) -> Iterator[AdditiveSequence[Element]]:
    """Yield every distinct atom divisor within the configured bound."""

    maximum_length = min(
        len(sequence),
        sequence.parent().davenport_bound,
    )
    zero = sequence.parent().base_parent.zero()
    for candidate in sequence.subsequences(
        nonempty=True,
        max_length=maximum_length,
    ):
        if candidate.is_atom() and zero not in candidate:
            yield candidate


@dataclass(frozen=True, slots=True)
class FactorizationStatistics:
    """Size statistics for one memoized factorization problem."""

    candidate_atoms: int
    states: int
    transitions: int


@dataclass(frozen=True, slots=True)
class _IndexedAtom(Generic[Element]):
    sequence: AdditiveSequence[Element]
    sparse_counts: tuple[tuple[int, int], ...]


class FactorizationSolver(Generic[Element]):
    """Solve one reduced factorization problem on immutable multiplicity states.

    Relevant nonzero atoms are discovered once, or filtered from a complete
    :class:`AtomCatalogue`, and encoded as sparse count vectors. The input
    must contain no identity terms; equal remainder sequences become one state
    in an acyclic graph.
    """

    def __init__(
        self,
        sequence: AdditiveSequence[Element],
        *,
        atom_catalogue: AtomCatalogue[Element] | None = None,
    ) -> None:
        if not isinstance(sequence, AdditiveSequence):
            raise TypeError("expected an additive sequence")
        if atom_catalogue is not None:
            if not isinstance(atom_catalogue, AtomCatalogue):
                raise TypeError("atom_catalogue must be an AtomCatalogue")
            if atom_catalogue.space is not sequence.parent():
                raise TypeError("atom catalogue and sequence use different spaces")
        if sequence.parent().base_parent.zero() in sequence:
            raise ValueError(
                "reduced factorizations require a sequence without zero terms"
            )

        self.sequence = sequence
        self.davenport_bound = sequence.parent().davenport_bound
        self.support = sequence.support
        self._support_index = {
            term: index for index, term in enumerate(self.support)
        }
        multiplicities = sequence.multiplicities
        self.initial_state: MultiplicityState = tuple(
            multiplicities[term] for term in self.support
        )
        self.zero_state: MultiplicityState = (0,) * len(self.support)

        candidates: Iterable[AdditiveSequence[Element]]
        if atom_catalogue is None:
            candidates = _atom_divisors(sequence)
        else:
            candidates = atom_catalogue.divisors(sequence)
        self._atoms = self._index_atoms(candidates)

        by_term: list[list[int]] = [[] for _ in self.support]
        for atom_index, atom in enumerate(self._atoms):
            for term_index, _ in atom.sparse_counts:
                by_term[term_index].append(atom_index)
        self._atoms_by_term = tuple(tuple(indices) for indices in by_term)
        self._transitions: dict[
            MultiplicityState, tuple[tuple[int, MultiplicityState], ...]
        ] | None = None
        self._length_bits: dict[MultiplicityState, int] | None = None
        self._fixed_length_results: dict[
            tuple[MultiplicityState, int], bool
        ] = {}
        factor_lengths = tuple(len(atom.sequence) for atom in self._atoms)
        self._minimum_factor_length = min(factor_lengths, default=0)
        self._maximum_factor_length = max(factor_lengths, default=0)
        self._minimum_length: int | None = None
        self._maximum_length: int | None = None
        self._minimum_length_is_solved = False
        self._maximum_length_is_solved = False

    def _index_atoms(
        self, atoms: Iterable[AdditiveSequence[Element]]
    ) -> tuple[_IndexedAtom[Element], ...]:
        indexed_by_counts: dict[MultiplicityState, _IndexedAtom[Element]] = {}
        for atom in atoms:
            if atom.parent() is not self.sequence.parent():
                raise TypeError("candidate atom and sequence use different spaces")
            counts = [0] * len(self.support)
            valid = True
            for term, count in atom.multiplicities.items():
                index = self._support_index.get(term)
                if index is None or count > self.initial_state[index]:
                    valid = False
                    break
                counts[index] = count
            if not valid or not atom:
                continue
            count_tuple = tuple(counts)
            indexed_by_counts[count_tuple] = _IndexedAtom(
                sequence=atom,
                sparse_counts=tuple(
                    (index, count)
                    for index, count in enumerate(count_tuple)
                    if count
                ),
            )
        return tuple(
            sorted(
                indexed_by_counts.values(),
                key=lambda atom: (len(atom.sequence), tuple(atom.sequence)),
            )
        )

    def _state_transitions(
        self, state: MultiplicityState
    ) -> tuple[tuple[int, MultiplicityState], ...]:
        if state == self.zero_state:
            return ()
        first_term = next(index for index, count in enumerate(state) if count)
        transitions = []
        for atom_index in self._atoms_by_term[first_term]:
            atom = self._atoms[atom_index]
            if all(state[index] >= count for index, count in atom.sparse_counts):
                remainder = list(state)
                for index, count in atom.sparse_counts:
                    remainder[index] -= count
                transitions.append((atom_index, tuple(remainder)))
        return tuple(transitions)

    def _query_transitions(
        self,
        state: MultiplicityState,
    ) -> tuple[tuple[int, MultiplicityState], ...]:
        if self._transitions is not None:
            known = self._transitions.get(state)
            if known is not None:
                return known
        return self._state_transitions(state)

    def _known_fixed_length_result(
        self,
        state: MultiplicityState,
        factor_count: int,
    ) -> bool | None:
        key = state, factor_count
        if key in self._fixed_length_results:
            return self._fixed_length_results[key]
        if self._length_bits is not None:
            bits = self._length_bits.get(state)
            if bits is not None:
                result = factor_count < bits.bit_length() and bool(
                    bits & (1 << factor_count)
                )
                self._fixed_length_results[key] = result
                return result
        if factor_count == 0:
            result = state == self.zero_state
        elif state == self.zero_state or not self._atoms:
            result = False
        else:
            terms = sum(state)
            result = None
            if (
                terms < self._minimum_factor_length * factor_count
                or terms > self._maximum_factor_length * factor_count
            ):
                result = False
        if result is not None:
            self._fixed_length_results[key] = result
        return result

    def _has_fixed_length_from_state(
        self,
        state: MultiplicityState,
        factor_count: int,
    ) -> bool:
        known = self._known_fixed_length_result(state, factor_count)
        if known is not None:
            return known

        frames: list[
            tuple[
                MultiplicityState,
                int,
                tuple[tuple[int, MultiplicityState], ...],
                int,
            ]
        ] = [(state, factor_count, self._query_transitions(state), 0)]
        while frames:
            current, remaining, transitions, next_transition = frames[-1]
            key = current, remaining
            if next_transition == len(transitions):
                self._fixed_length_results[key] = False
                frames.pop()
                continue

            frames[-1] = (
                current,
                remaining,
                transitions,
                next_transition + 1,
            )
            _, remainder = transitions[next_transition]
            next_remaining = remaining - 1
            child_result = self._known_fixed_length_result(
                remainder,
                next_remaining,
            )
            if child_result is True:
                self._fixed_length_results[key] = True
                for ancestor, ancestor_remaining, _, _ in frames:
                    self._fixed_length_results[
                        ancestor, ancestor_remaining
                    ] = True
                return True
            if child_result is False:
                continue
            frames.append(
                (
                    remainder,
                    next_remaining,
                    self._query_transitions(remainder),
                    0,
                )
            )
        return False

    @staticmethod
    def _extreme_lengths_from_bits(bits: int) -> tuple[int | None, int | None]:
        if not bits:
            return None, None
        minimum = (bits & -bits).bit_length() - 1
        maximum = bits.bit_length() - 1
        return minimum, maximum

    def _build_state_graph(
        self,
    ) -> dict[MultiplicityState, tuple[tuple[int, MultiplicityState], ...]]:
        if self._transitions is None:
            transitions: dict[
                MultiplicityState, tuple[tuple[int, MultiplicityState], ...]
            ] = {}
            pending = [self.initial_state]
            discovered = {self.initial_state}
            while pending:
                state = pending.pop()
                state_transitions = self._state_transitions(state)
                transitions[state] = state_transitions
                for _, remainder in state_transitions:
                    if remainder not in discovered:
                        discovered.add(remainder)
                        pending.append(remainder)
            self._transitions = transitions
        return self._transitions

    def _solve_lengths(self) -> dict[MultiplicityState, int]:
        if self._length_bits is None:
            transitions = self._build_state_graph()
            length_bits: dict[MultiplicityState, int] = {}
            for state in sorted(transitions, key=sum):
                if state == self.zero_state:
                    length_bits[state] = 1
                    continue
                bits = 0
                for _, remainder in transitions[state]:
                    bits |= length_bits[remainder] << 1
                length_bits[state] = bits
            self._length_bits = length_bits
        return self._length_bits

    def _solve_extreme_lengths(self) -> tuple[int | None, int | None]:
        if self._minimum_length_is_solved and self._maximum_length_is_solved:
            return self._minimum_length, self._maximum_length

        if self._length_bits is not None:
            minimum, maximum = self._extreme_lengths_from_bits(
                self._length_bits[self.initial_state]
            )
        else:
            transitions = self._build_state_graph()
            minimum_by_state: dict[MultiplicityState, int | None] = {}
            maximum_by_state: dict[MultiplicityState, int | None] = {}
            for state in sorted(transitions, key=sum):
                if state == self.zero_state:
                    minimum_by_state[state] = 0
                    maximum_by_state[state] = 0
                    continue
                child_bounds = [
                    (minimum_by_state[remainder], maximum_by_state[remainder])
                    for _, remainder in transitions[state]
                    if minimum_by_state[remainder] is not None
                ]
                if child_bounds:
                    minimum_by_state[state] = 1 + min(
                        minimum
                        for minimum, _ in child_bounds
                        if minimum is not None
                    )
                    maximum_by_state[state] = 1 + max(
                        maximum
                        for _, maximum in child_bounds
                        if maximum is not None
                    )
                else:
                    minimum_by_state[state] = None
                    maximum_by_state[state] = None
            minimum = minimum_by_state[self.initial_state]
            maximum = maximum_by_state[self.initial_state]

        self._minimum_length = minimum
        self._maximum_length = maximum
        self._minimum_length_is_solved = True
        self._maximum_length_is_solved = True
        return minimum, maximum

    def _sequence_from_state(
        self, state: MultiplicityState
    ) -> AdditiveSequence[Element]:
        return self.sequence.parent().from_multiplicities(
            {
                term: state[index]
                for index, term in enumerate(self.support)
                if state[index]
            }
        )

    @property
    def statistics(self) -> FactorizationStatistics:
        """Return atom, state and transition counts for this problem."""

        transitions = self._build_state_graph()
        return FactorizationStatistics(
            candidate_atoms=len(self._atoms),
            states=len(transitions),
            transitions=sum(len(edges) for edges in transitions.values()),
        )

    def length_set(self) -> set[int]:
        """Return the complete set of attained factorization lengths."""

        length_bits = self._solve_lengths()
        bits = length_bits[self.initial_state]
        return {
            length
            for length in range(bits.bit_length())
            if bits & (1 << length)
        }

    def minimum_factorization_length(self) -> int | None:
        """Return the minimum factorization length, or ``None`` if absent."""

        if self._minimum_length_is_solved:
            return self._minimum_length
        if self._length_bits is not None:
            minimum, _ = self._extreme_lengths_from_bits(
                self._length_bits[self.initial_state]
            )
            self._minimum_length = minimum
            self._minimum_length_is_solved = True
            return minimum

        pending = deque(((self.initial_state, 0),))
        discovered = {self.initial_state}
        while pending:
            state, length = pending.popleft()
            if state == self.zero_state:
                self._minimum_length = length
                self._minimum_length_is_solved = True
                return length
            for _, remainder in self._query_transitions(state):
                if remainder not in discovered:
                    discovered.add(remainder)
                    pending.append((remainder, length + 1))

        self._minimum_length = None
        self._minimum_length_is_solved = True
        return None

    def maximum_factorization_length(self) -> int | None:
        """Return the maximum factorization length, or ``None`` if absent."""

        if self._maximum_length_is_solved:
            return self._maximum_length
        _, maximum = self._solve_extreme_lengths()
        return maximum

    def has_factorization_of_length(self, factor_count: int) -> bool:
        """Return whether a factorization has exactly ``factor_count`` factors."""

        factor_count = _non_negative_integer(factor_count, name="factor count")
        return self._has_fixed_length_from_state(
            self.initial_state,
            factor_count,
        )

    def factorization_witnesses(self) -> dict[int, Factorization]:
        """Return one deterministic factorization per attained length."""

        length_bits = self._solve_lengths()
        transitions = self._build_state_graph()
        witnesses = {}
        for length in sorted(self.length_set()):
            state = self.initial_state
            remaining_length = length
            factors = []
            while remaining_length:
                for atom_index, remainder in transitions[state]:
                    if length_bits[remainder] & (1 << (remaining_length - 1)):
                        factors.append(self._atoms[atom_index].sequence)
                        state = remainder
                        remaining_length -= 1
                        break
                else:  # pragma: no cover - guards the dynamic-programming invariant
                    raise RuntimeError("could not reconstruct factorization witness")
            witnesses[length] = tuple(factors)
        return witnesses

    def factorization_witness(
        self,
        factor_count: int,
        *,
        minimum_matching_factors: int = 0,
        factor_predicate: Callable[[AdditiveSequence[Element]], bool] | None = None,
    ) -> Factorization | None:
        """Find a factorization of a requested length and optional profile.

        When ``minimum_matching_factors`` is positive, at least that many
        factors must satisfy ``factor_predicate``. The search reuses the
        remainder DAG and its length bitsets rather than enumerating every
        factorization.
        """

        factor_count = _non_negative_integer(factor_count, name="factor count")
        minimum_matching_factors = _non_negative_integer(
            minimum_matching_factors,
            name="minimum matching factors",
        )
        if minimum_matching_factors > factor_count:
            raise ValueError(
                "minimum matching factors cannot exceed factor count"
            )
        if minimum_matching_factors and not callable(factor_predicate):
            raise ValueError(
                "factor_predicate is required when matching factors are required"
            )

        length_bits = self._solve_lengths()
        transitions = self._build_state_graph()
        if not length_bits[self.initial_state] & (1 << factor_count):
            return None

        frames: list[tuple[MultiplicityState, int, int, int]] = [
            (self.initial_state, factor_count, minimum_matching_factors, 0)
        ]
        path: list[int] = []
        dead: set[tuple[MultiplicityState, int, int]] = set()

        while frames:
            state, remaining, required_matches, next_transition = frames[-1]
            key = (state, remaining, required_matches)

            if remaining == 0:
                if state == self.zero_state and required_matches == 0:
                    return tuple(self._atoms[index].sequence for index in path)
                dead.add(key)
                frames.pop()
                if path:
                    path.pop()
                continue

            state_transitions = transitions[state]
            if next_transition >= len(state_transitions):
                dead.add(key)
                frames.pop()
                if path:
                    path.pop()
                continue

            frames[-1] = (
                state,
                remaining,
                required_matches,
                next_transition + 1,
            )
            atom_index, remainder = state_transitions[next_transition]
            atom = self._atoms[atom_index].sequence
            matches = factor_predicate is not None and factor_predicate(atom)
            next_required_matches = max(0, required_matches - int(matches))
            next_remaining = remaining - 1
            next_key = (remainder, next_remaining, next_required_matches)
            if next_required_matches > next_remaining or next_key in dead:
                continue
            if not length_bits[remainder] & (1 << next_remaining):
                dead.add(next_key)
                continue

            path.append(atom_index)
            frames.append(
                (remainder, next_remaining, next_required_matches, 0)
            )

        return None

    def factorizations(
        self,
        *,
        factor_count: int | None = None,
    ) -> Iterator[Factorization]:
        """Yield each unordered factorization, optionally of one length."""

        if factor_count is not None:
            factor_count = _non_negative_integer(factor_count, name="factor count")
            if not self.has_factorization_of_length(factor_count):
                return

        path: list[AdditiveSequence[Element]] = []
        stack: list[tuple[MultiplicityState, int]] = [(self.initial_state, 0)]
        while stack:
            state, next_atom = stack[-1]
            if state == self.zero_state:
                if factor_count is None or len(path) == factor_count:
                    yield tuple(path)
                stack.pop()
                if stack:
                    path.pop()
                continue

            if factor_count is not None:
                remaining_factors = factor_count - len(path)
                if not self._has_fixed_length_from_state(
                    state,
                    remaining_factors,
                ):
                    stack.pop()
                    if stack:
                        path.pop()
                    continue

            if next_atom == len(self._atoms):
                stack.pop()
                if stack:
                    path.pop()
                continue

            stack[-1] = (state, next_atom + 1)
            atom = self._atoms[next_atom]
            if all(state[index] >= count for index, count in atom.sparse_counts):
                remainder = list(state)
                for index, count in atom.sparse_counts:
                    remainder[index] -= count
                path.append(atom.sequence)
                stack.append((tuple(remainder), next_atom))

    def digraph(self):
        """Return the memoized remainder graph as a NetworkX ``DiGraph``.

        Each edge stores the removed atom in its ``"atom"`` attribute.
        """

        transitions_by_state = self._build_state_graph()
        sequences = {
            state: self._sequence_from_state(state)
            for state in transitions_by_state
        }
        graph = nx.DiGraph()
        graph.add_nodes_from(sequences.values())
        for state, transitions in transitions_by_state.items():
            for atom_index, remainder in transitions:
                graph.add_edge(
                    sequences[state],
                    sequences[remainder],
                    atom=self._atoms[atom_index].sequence,
                )
        return graph
