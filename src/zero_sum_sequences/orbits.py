"""Automorphism actions on additive sequences.

The core orbit routines are deliberately independent of Sage or any other
computer algebra system.  An :class:`AutomorphismAction` consists of a finite
list of generators and, optionally, a callback describing how a generator
acts on one base-parent element.  For the common case in which generators are
callables, the callback is inferred automatically.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass, field
from operator import index
from typing import Generic, TypeVar

from .additive_sequence import AdditiveSequence

Element = TypeVar("Element")
Generator = TypeVar("Generator")


class AutomorphismActionUnavailable(RuntimeError):
    """Raised when an additive parent cannot supply an automorphism action."""


def _non_negative_index(value: object, *, name: str) -> int:
    """Return a non-negative indexable integer, rejecting booleans."""

    if isinstance(value, bool):
        raise ValueError(f"{name} must be a non-negative integer")
    try:
        integer = index(value)
    except TypeError:
        raise ValueError(f"{name} must be a non-negative integer") from None
    if integer < 0:
        raise ValueError(f"{name} must be a non-negative integer")
    return integer


class AutomorphismAction(Generic[Element, Generator]):
    """A finite generating set acting on terms of an additive parent.

    Parameters
    ----------
    generators:
        The generators, in the order used by the orbit search.  If no
        ``apply_term`` callback is supplied, every generator must be callable
        on one base-parent element.
    apply_term:
        Optional callback ``apply_term(generator, term)``.  This is useful for
        matrix or CAS objects that represent maps without being callable
        themselves.
    materialize_word:
        Optional callback turning an :class:`OrbitWitness` into a concrete
        automorphism.  Orbit traversal itself only needs generator words;
        adapters can use this hook when an explicit map or matrix is needed.

    The class intentionally does not inspect or enumerate the generated group.
    It only provides the term-wise action needed by the orbit routines.
    """

    __slots__ = ("generators", "_apply_term", "_materialize_word")

    def __init__(
        self,
        generators: Iterable[Generator],
        *,
        apply_term: Callable[[Generator, Element], Element] | None = None,
        materialize_word: Callable[[OrbitWitness], object] | None = None,
    ) -> None:
        try:
            normalized_generators = tuple(generators)
        except TypeError:
            raise TypeError("automorphism generators must be iterable") from None
        if apply_term is not None and not callable(apply_term):
            raise TypeError("apply_term must be callable")
        if materialize_word is not None and not callable(materialize_word):
            raise TypeError("materialize_word must be callable")
        if apply_term is None:
            for generator in normalized_generators:
                if not callable(generator):
                    raise TypeError(
                        "generators must be callable when apply_term is omitted"
                    )
            apply_term = _call_generator
        self.generators = normalized_generators
        self._apply_term = apply_term
        self._materialize_word = materialize_word

    def __iter__(self) -> Iterator[Generator]:
        return iter(self.generators)

    def __len__(self) -> int:
        return len(self.generators)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.generators!r})"

    def apply_term(self, generator_index: int, term: Element) -> Element:
        """Apply one indexed generator to a base-parent term."""

        generator_index = _non_negative_index(
            generator_index,
            name="generator index",
        )
        try:
            generator = self.generators[generator_index]
        except IndexError:
            raise IndexError("automorphism generator index out of range") from None
        return self._apply_term(generator, term)

    def apply_sequence(
        self,
        sequence: AdditiveSequence[Element],
        generator_index: int,
    ) -> AdditiveSequence[Element]:
        """Apply one indexed generator term by term to ``sequence``."""

        if not isinstance(sequence, AdditiveSequence):
            raise TypeError("expected an additive sequence")
        return sequence.map_terms(
            lambda term: self.apply_term(generator_index, term)
        )

    def apply_word(
        self,
        sequence: AdditiveSequence[Element],
        word: OrbitWitness | Iterable[int],
    ) -> AdditiveSequence[Element]:
        """Apply a generator word from left to right to ``sequence``."""

        if not isinstance(word, OrbitWitness):
            word = OrbitWitness(word)
        image = sequence
        for generator_index in word.generator_indices:
            image = self.apply_sequence(image, generator_index)
        return image

    def materialize_word(
        self,
        word: OrbitWitness | Iterable[int],
    ):
        """Return the concrete automorphism represented by ``word``.

        Not every action has a useful concrete representation.  Such actions
        still support orbit traversal, but this method raises ``TypeError``.
        """

        if not isinstance(word, OrbitWitness):
            word = OrbitWitness(word)
        if self._materialize_word is None:
            raise TypeError("this action cannot materialize generator words")
        return self._materialize_word(word)


def _call_generator(generator, term):
    return generator(term)


@dataclass(frozen=True, slots=True, init=False)
class OrbitWitness:
    """A word in the indexed generators of an :class:`AutomorphismAction`.

    ``generator_indices`` records the order in which generators are applied.
    The empty tuple is the identity word.
    """

    generator_indices: tuple[int, ...]
    _sequence_space: object | None = field(
        default=None,
        repr=False,
        compare=False,
        hash=False,
    )
    _action: AutomorphismAction | None = field(
        default=None,
        repr=False,
        compare=False,
        hash=False,
    )

    def __init__(
        self,
        generator_indices: Iterable[int] = (),
        *,
        _sequence_space=None,
        _action: AutomorphismAction | None = None,
    ) -> None:
        try:
            normalized = tuple(
                _non_negative_index(value, name="generator indices")
                for value in generator_indices
            )
        except TypeError:
            raise TypeError("generator indices must be iterable") from None
        object.__setattr__(self, "generator_indices", normalized)
        object.__setattr__(self, "_sequence_space", _sequence_space)
        object.__setattr__(self, "_action", _action)

    def __iter__(self) -> Iterator[int]:
        return iter(self.generator_indices)

    def __len__(self) -> int:
        return len(self.generator_indices)

    def __repr__(self) -> str:
        return (
            f"{type(self).__name__}("
            f"generator_indices={self.generator_indices!r})"
        )

    def _with_context(self, sequence_space, action) -> OrbitWitness:
        return type(self)(
            self.generator_indices,
            _sequence_space=sequence_space,
            _action=action,
        )

    def additive_generator_images(self) -> tuple[tuple[object, object], ...]:
        """Return the induced images of distinguished additive generators.

        Witnesses returned by :func:`orbit_witness` remember the sequence
        space and action used to construct them.  The base parent must expose
        an additive generating set through ``additive_generators()``,
        ``basis()``, or ``gens()``.
        """

        if self._sequence_space is None or self._action is None:
            raise ValueError(
                "term images require a witness returned by orbit_witness()"
            )
        parent = self._sequence_space.base_parent
        generators = _additive_generators(parent)
        images = []
        for generator in generators:
            image = generator
            for generator_index in self.generator_indices:
                image = self._action.apply_term(generator_index, image)
            images.append((generator, image))
        return tuple(images)

    def show(self, *, file=None) -> None:
        """Print the induced map on distinguished additive generators."""

        for generator, image in self.additive_generator_images():
            print(f"{generator} ↦ {image}", file=file)


def _additive_generators(parent) -> tuple:
    for provider_name in ("additive_generators", "basis", "gens"):
        provider = getattr(parent, provider_name, None)
        if not callable(provider):
            continue
        try:
            return tuple(provider())
        except (AttributeError, NotImplementedError, TypeError, ValueError):
            continue
    raise ValueError(
        "the additive parent does not expose a basis or generating set"
    )


def _action_from_parent_provider(parent) -> AutomorphismAction | None:
    """Return an action explicitly exposed by ``parent``, if available."""

    provider = getattr(parent, "automorphism_action", None)
    if isinstance(provider, AutomorphismAction):
        return provider
    if callable(provider):
        try:
            action = provider()
        except (AttributeError, NotImplementedError):
            action = None
        if action is not None:
            if not isinstance(action, AutomorphismAction):
                raise TypeError(
                    "automorphism_action() must return an AutomorphismAction"
                )
            return action

    provider = getattr(parent, "automorphism_generators", None)
    if provider is not None:
        try:
            generators = provider() if callable(provider) else provider
        except (AttributeError, NotImplementedError):
            generators = None
        if generators is not None:
            return AutomorphismAction(generators)
    return None


def _sage_vector_space_action(parent) -> AutomorphismAction | None:
    """Construct the natural matrix action for a finite Sage vector space."""

    try:
        from sage.all import FiniteFields, GL, VectorSpaces, identity_matrix, matrix
    except ImportError:
        return None

    try:
        field = parent.base_ring()
        if field not in FiniteFields():
            return None
        if parent not in VectorSpaces(field).FiniteDimensional():
            return None
        dimension = int(parent.dimension())
    except (AttributeError, TypeError, ValueError):
        return None

    generators = []
    for generator in GL(dimension, field).gens() if dimension else ():
        generator = matrix(generator)
        generator.set_immutable()
        generators.append(generator)
    generators = tuple(generators)

    def apply_term(generator, term):
        return parent(generator * parent(term))

    def materialize_word(word):
        result = identity_matrix(field, dimension)
        for generator_index in word.generator_indices:
            # Generator words record maps in application order.  With column
            # vectors, the matrix of g_0 followed by g_1 is g_1 * g_0.
            result = generators[generator_index] * result
        result = matrix(result)
        result.set_immutable()
        return result

    return AutomorphismAction(
        generators,
        apply_term=apply_term,
        materialize_word=materialize_word,
    )


def _cas_automorphism_action(parent) -> AutomorphismAction | None:
    """Try conservative automorphism-group conventions used by CAS parents."""

    is_finite = getattr(parent, "is_finite", None)
    try:
        if not callable(is_finite) or not is_finite():
            return None
    except (NotImplementedError, TypeError, ValueError):
        return None

    for method_name in ("automorphism_group", "aut"):
        method = getattr(parent, method_name, None)
        if not callable(method):
            continue
        try:
            group = method()
            generators = tuple(group.gens())
        except (AttributeError, NotImplementedError, TypeError, ValueError):
            continue
        if not all(callable(generator) for generator in generators):
            continue

        def apply_term(generator, term):
            return parent(generator(term))

        action = AutomorphismAction(generators, apply_term=apply_term)
        try:
            if any(
                action.apply_term(index, parent.zero()) != parent.zero()
                for index in range(len(action))
            ):
                continue
        except (TypeError, ValueError):
            continue
        return action
    return None


def automorphism_action(sequence_space) -> AutomorphismAction:
    """Resolve and cache an automorphism action for a sequence space.

    Resolution first honors an explicit provider on the additive parent, then
    recognizes finite-dimensional Sage vector spaces over finite fields, and
    finally tries conventional CAS automorphism-group methods.
    """

    cached = sequence_space._automorphism_action
    if cached is not None:
        return cached

    parent = sequence_space.base_parent
    action = _action_from_parent_provider(parent)
    if action is None:
        action = _sage_vector_space_action(parent)
    if action is None:
        action = _cas_automorphism_action(parent)
    if action is None:
        raise AutomorphismActionUnavailable(
            "could not determine automorphism generators for the additive "
            "parent; pass action=AutomorphismAction(...) explicitly"
        )
    sequence_space._automorphism_action = action
    return action


def _normalize_action(action, sequence) -> AutomorphismAction:
    if isinstance(action, AutomorphismAction):
        return action
    if action is None:
        return automorphism_action(sequence.parent())
    return AutomorphismAction(action)


def _require_sequence_pair(
    source: AdditiveSequence[Element],
    target: AdditiveSequence[Element],
) -> None:
    if not isinstance(source, AdditiveSequence):
        raise TypeError("expected an additive sequence")
    source._require_same_space(target)


def _orbit_with_witnesses(
    source: AdditiveSequence[Element],
    action: AutomorphismAction,
    *,
    target: AdditiveSequence[Element] | None = None,
) -> tuple[dict[AdditiveSequence[Element], OrbitWitness], bool]:
    """Traverse an orbit and optionally stop when ``target`` is found."""

    witnesses: dict[AdditiveSequence[Element], OrbitWitness] = {
        source: OrbitWitness()
    }
    pending = deque((source,))
    while pending:
        current = pending.popleft()
        if target is not None and current == target:
            return witnesses, True
        current_word = witnesses[current]
        for generator_index in range(len(action)):
            image = action.apply_sequence(current, generator_index)
            if image in witnesses:
                continue
            image_word = OrbitWitness(
                (*current_word.generator_indices, generator_index)
            )
            witnesses[image] = image_word
            if target is not None and image == target:
                return witnesses, True
            pending.append(image)
    return witnesses, target is None


def orbit(
    sequence: AdditiveSequence[Element],
    *,
    action: AutomorphismAction[Element, Generator] | Iterable[Generator] | None = None,
) -> tuple[AdditiveSequence[Element], ...]:
    """Return all distinct automorphism images of ``sequence``.

    Images are returned in deterministic breadth-first order, beginning with
    ``sequence`` itself.  The action must have a finite orbit for this
    materialized traversal to terminate.
    """

    if not isinstance(sequence, AdditiveSequence):
        raise TypeError("expected an additive sequence")
    normalized_action = _normalize_action(action, sequence)
    witnesses, _ = _orbit_with_witnesses(sequence, normalized_action)
    return tuple(witnesses)


def is_in_same_orbit(
    source: AdditiveSequence[Element],
    target: AdditiveSequence[Element],
    *,
    action: AutomorphismAction[Element, Generator] | Iterable[Generator] | None = None,
) -> bool:
    """Return whether two sequences lie in one automorphism orbit.

    When ``action`` is omitted, equality is recognized without resolving an
    action from the additive parent.
    """

    _require_sequence_pair(source, target)
    if source == target and action is None:
        return True
    normalized_action = _normalize_action(action, source)
    _, found = _orbit_with_witnesses(
        source,
        normalized_action,
        target=target,
    )
    return found


def orbit_witness(
    source: AdditiveSequence[Element],
    target: AdditiveSequence[Element],
    *,
    action: AutomorphismAction[Element, Generator] | Iterable[Generator] | None = None,
) -> OrbitWitness | None:
    """Return a shortest deterministic generator word to ``target``.

    The empty word is returned for equal sequences.  If ``action`` is omitted,
    this equality case does not require the additive parent to expose an
    action.  If ``target`` is not in the orbit, return ``None``.
    """

    _require_sequence_pair(source, target)
    if source == target and action is None:
        return OrbitWitness()._with_context(
            source.parent(),
            AutomorphismAction(()),
        )
    normalized_action = _normalize_action(action, source)
    witnesses, found = _orbit_with_witnesses(
        source,
        normalized_action,
        target=target,
    )
    if not found:
        return None
    return witnesses[target]._with_context(source.parent(), normalized_action)


__all__ = [
    "AutomorphismAction",
    "AutomorphismActionUnavailable",
    "OrbitWitness",
    "automorphism_action",
    "is_in_same_orbit",
    "orbit",
    "orbit_witness",
]
