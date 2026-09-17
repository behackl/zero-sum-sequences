"""Materialized automorphism groups, stabilizers and canonical forms.

An :class:`AutomorphismGroup` is a finite group of automorphisms of an
additive parent, given by *explicit elements* together with callbacks that
apply, compose and invert them.  The parent delegates: it may construct the
group itself (:meth:`FiniteAdditiveGroup.automorphism_group`), a CAS adapter
may supply matrices, or the group is obtained as the closure of a finite
generating set of term maps.

The generic element type is :class:`Automorphism`, an automorphism of a
finite parent stored as its *image table*: for every parent element (in the
parent's canonical order) the index of its image.  It needs nothing but a
finite iterable parent, so it works for every :class:`FiniteAdditiveGroup`.

The group acts on additive sequences and on tuples of them.  Tuples are
treated as multisets unless ``ordered=True`` is requested, which is the
natural notion for contexts, decompositions and unordered pairs.  Canonical
forms are minima of orbits with respect to the total order on sequences
(length first, then the sorted term list), and the automorphism returned
with a canonical form is the smallest group element realising it, so that
witnesses do not depend on the order in which the group was discovered.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass
from functools import partial
from typing import Generic, TypeVar

from .additive_sequence import AdditiveSequence, AdditiveSequenceSpace, _immutable_term

Element = TypeVar("Element")
Auto = TypeVar("Auto")


class AutomorphismGroupUnavailable(RuntimeError):
    """Raised when no materialized automorphism group can be obtained."""


# Image-table automorphisms of a finite parent


class _ParentTable(Generic[Element]):
    """Shared lookup data for the image-table automorphisms of one parent."""

    __slots__ = ("elements", "index", "generator_indices")

    def __init__(self, parent) -> None:
        try:
            elements = tuple(_immutable_term(term) for term in parent)
        except TypeError:
            raise TypeError(
                "image-table automorphisms require a finite iterable parent"
            ) from None
        try:
            self.elements = tuple(sorted(elements))
        except TypeError:
            self.elements = elements
        self.index = {term: position for position, term in enumerate(self.elements)}
        if len(self.index) != len(self.elements):
            raise ValueError("the parent contains repeated elements")
        from .orbits import _additive_generators

        try:
            generators = tuple(_immutable_term(g) for g in _additive_generators(parent))
        except (NotImplementedError, TypeError, ValueError):
            generators = ()
        self.generator_indices = tuple(
            self.index[g] for g in generators if g in self.index
        )


class Automorphism(Generic[Element]):
    """An automorphism of a finite parent stored as an image table.

    Instances compare equal, hash and order by their image table, which makes
    the element order of a group canonical.  ``generator_images`` gives the
    compact description by the images of the parent's additive generators
    (for a coordinate group this is the matrix, column by column).
    """

    __slots__ = ("_table", "_data", "_hash")

    def __init__(self, table: tuple[int, ...], data: _ParentTable[Element]) -> None:
        self._table = table
        self._data = data
        self._hash = hash(table)

    @property
    def table(self) -> tuple[int, ...]:
        """Image indices of the parent elements, in parent order."""

        return self._table

    def __call__(self, term: Element) -> Element:
        data = self._data
        return data.elements[self._table[data.index[_immutable_term(term)]]]

    def compose(self, other: Automorphism[Element]) -> Automorphism[Element]:
        """Return ``self ∘ other`` (apply ``other`` first)."""

        table = self._table
        return Automorphism(tuple(table[i] for i in other._table), self._data)

    def inverse(self) -> Automorphism[Element]:
        inverse = [0] * len(self._table)
        for position, image in enumerate(self._table):
            inverse[image] = position
        return Automorphism(tuple(inverse), self._data)

    def is_identity(self) -> bool:
        return all(position == image for position, image in enumerate(self._table))

    def generator_images(self) -> tuple[Element, ...]:
        """Images of the parent's additive generators, if it exposes any."""

        data = self._data
        return tuple(data.elements[self._table[i]] for i in data.generator_indices)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Automorphism) and self._table == other._table

    def __hash__(self) -> int:
        return self._hash

    def __lt__(self, other: Automorphism[Element]) -> bool:
        return self._table < other._table

    def __repr__(self) -> str:
        images = self.generator_images()
        if images:
            return f"Automorphism(generator_images={images!r})"
        return f"Automorphism(table={self._table!r})"


def _table_of(term_map: Callable[[Element], Element], data: _ParentTable) -> tuple[int, ...]:
    try:
        table = tuple(data.index[_immutable_term(term_map(term))] for term in data.elements)
    except KeyError:
        raise ValueError("a generator maps a parent element outside the parent") from None
    if len(set(table)) != len(table):
        raise ValueError("a generator is not a bijection of the parent")
    return table


# The group


@dataclass(frozen=True)
class CanonicalForm(Generic[Auto]):
    """Result of :meth:`AutomorphismGroup.canonical_form`.

    ``image`` is the orbit minimum, ``automorphism`` the smallest group
    element with ``apply(automorphism, obj) == image``, and for tuple
    objects ``permutation`` records where each input position ends up in
    the (sorted) image: ``image[permutation[i]] == apply(automorphism, obj[i])``.
    """

    image: object
    automorphism: Auto
    permutation: tuple[int, ...] | None = None


class AutomorphismGroup(Generic[Element, Auto]):
    """A finite group of automorphisms given by explicit elements.

    Parameters
    ----------
    elements:
        The group elements.  They must be hashable; they are stored sorted by
        ``key`` (default: the elements' own order), identity first.
    apply_term:
        ``apply_term(auto, term)``: the action on one parent element.
    compose, inverse, identity:
        Group structure.  ``compose(a, b)`` is ``a ∘ b``.  Only needed for
        :meth:`compose`/:meth:`inverse` and for word witnesses; orbits,
        stabilizers and canonical forms use ``apply_term`` alone.
    generators:
        A generating subset used for orbit traversals, whose cost is then
        proportional to the orbit size.  Defaults to all elements.
    """

    def __init__(
        self,
        elements: Iterable[Auto],
        *,
        apply_term: Callable[[Auto, Element], Element],
        compose: Callable[[Auto, Auto], Auto] | None = None,
        inverse: Callable[[Auto], Auto] | None = None,
        identity: Auto | None = None,
        generators: Iterable[Auto] | None = None,
        key: Callable[[Auto], object] | None = None,
    ) -> None:
        if not callable(apply_term):
            raise TypeError("apply_term must be callable")
        unique = list(dict.fromkeys(elements))
        if not unique:
            raise ValueError("an automorphism group needs at least the identity")
        self._key = key
        ordered = tuple(sorted(unique, key=key))
        self._element_set = frozenset(ordered)
        self._apply_term = apply_term
        self._compose = compose
        self._inverse = inverse
        if identity is None:
            identity = next((e for e in ordered if _is_identity(e)), None)
            if identity is None:
                raise ValueError("the elements do not contain the identity")
        elif identity not in self._element_set:
            raise ValueError("the identity is not among the elements")
        self._identity = identity
        # identity first, then the canonical order
        self._elements: tuple[Auto, ...] = (
            identity,
            *(e for e in ordered if e != identity),
        )
        if generators is None:
            self._generators = tuple(e for e in self._elements if e != identity)
        else:
            self._generators = tuple(dict.fromkeys(generators))
            if not set(self._generators) <= self._element_set:
                raise ValueError("generators must be elements of the group")
            if compose is not None and len(self._closure(self._generators)) != len(self):
                raise ValueError("the supplied generators do not generate the group")

    def _closure(self, generators: Iterable[Auto]) -> set[Auto]:
        """The subgroup generated by ``generators`` (needs ``compose``)."""

        seen = {self._identity}
        pending = [self._identity]
        while pending:
            current = pending.pop()
            for generator in generators:
                image = self._compose(generator, current)
                if image not in seen:
                    seen.add(image)
                    pending.append(image)
        return seen

    # structure

    @property
    def identity(self) -> Auto:
        return self._identity

    @property
    def elements(self) -> tuple[Auto, ...]:
        return self._elements

    @property
    def generators(self) -> tuple[Auto, ...]:
        return self._generators

    def __len__(self) -> int:
        return len(self._elements)

    def __iter__(self) -> Iterator[Auto]:
        return iter(self._elements)

    def __contains__(self, element: object) -> bool:
        return element in self._element_set

    def __repr__(self) -> str:
        return f"AutomorphismGroup(order={len(self)})"

    def compose(self, left: Auto, right: Auto) -> Auto:
        if self._compose is None:
            raise TypeError("this automorphism group cannot compose elements")
        return self._compose(left, right)

    def inverse(self, element: Auto) -> Auto:
        if self._inverse is None:
            raise TypeError("this automorphism group cannot invert elements")
        return self._inverse(element)

    def subgroup(self, elements: Iterable[Auto]) -> AutomorphismGroup[Element, Auto]:
        """A group on a subset of the elements, sharing all callbacks.

        When the group can compose, closure is checked and a small
        generating set of the subgroup is found greedily, so that orbit
        traversals in the subgroup stay cheap; otherwise the subset is not
        checked and all its non-identity elements serve as generators.
        """

        members = tuple(dict.fromkeys(elements))
        generators = None
        if self._compose is not None:
            generators = []
            closure = {self._identity}
            for element in sorted(members, key=self._key):
                if element not in closure:
                    generators.append(element)
                    closure = self._closure(generators)
            if len(closure) != len(members):
                raise ValueError("the elements are not closed under composition")
        return AutomorphismGroup(
            members,
            apply_term=self._apply_term,
            compose=self._compose,
            inverse=self._inverse,
            identity=self._identity,
            generators=generators,
            key=self._key,
        )

    # actions

    def apply_term(self, element: Auto, term: Element) -> Element:
        return self._apply_term(element, term)

    def apply(self, element: Auto, obj, *, ordered: bool = False):
        """Apply ``element`` to a sequence or to a tuple of sequences.

        Tuple images are sorted unless ``ordered``; ``apply_term`` must
        return elements of the parent.
        """

        if isinstance(obj, AdditiveSequence):
            return obj._map_canonical(
                lambda term: _immutable_term(self._apply_term(element, term))
            )
        if isinstance(obj, tuple):
            images = tuple(self.apply(element, item) for item in obj)
            return images if ordered else tuple(sorted(images))
        raise TypeError("expected an additive sequence or a tuple of them")

    def orbit_words(self, obj, *, ordered: bool = False) -> dict:
        """The orbit of ``obj`` as ``image -> generator word``.

        Words are tuples of indices into :attr:`generators`, applied left to
        right; ``obj`` itself has the empty word.  The traversal is
        breadth-first and costs proportionally to the orbit size.
        """

        obj = _prepare(obj, ordered)
        words = {obj: ()}
        pending = deque([obj])
        while pending:
            current = pending.popleft()
            word = words[current]
            for position, generator in enumerate(self._generators):
                image = self.apply(generator, current, ordered=ordered)
                if image not in words:
                    words[image] = (*word, position)
                    pending.append(image)
        return words

    def materialize(self, word: Iterable[int]) -> Auto:
        """The product of the generators of ``word`` (needs ``compose``)."""

        element = self._identity
        for position in word:
            element = self.compose(self._generators[position], element)
        return element

    def orbit(self, obj, *, ordered: bool = False) -> tuple:
        """All distinct images of ``obj``, sorted."""

        return tuple(sorted(self.orbit_words(obj, ordered=ordered)))

    def stabilizer(self, obj, *, ordered: bool = False) -> AutomorphismGroup[Element, Auto]:
        """The subgroup fixing ``obj`` (as a multiset for tuples)."""

        obj = _prepare(obj, ordered)
        return self.subgroup(
            e for e in self._elements if self.apply(e, obj, ordered=ordered) == obj
        )

    def transporter(self, source, target, *, ordered: bool = False) -> Auto | None:
        """The smallest element mapping ``source`` to ``target``, or ``None``."""

        source = _prepare(source, ordered)
        target = _prepare(target, ordered)
        for element in self._elements:
            if self.apply(element, source, ordered=ordered) == target:
                return element
        return None

    def canonical(self, obj, *, ordered: bool = False):
        """The orbit minimum of ``obj``."""

        return min(self.orbit_words(obj, ordered=ordered))

    def canonical_form(
        self, obj, *, ordered: bool = False, witness: str = "smallest"
    ) -> CanonicalForm[Auto]:
        """The orbit minimum together with an automorphism realising it.

        ``witness="smallest"`` scans the group for the smallest element
        mapping ``obj`` to the minimum (cost: the group order).
        ``witness="word"`` materializes the breadth-first generator word
        instead (cost: the orbit size; needs ``compose``); it is
        deterministic for a fixed generator list.
        """

        if witness not in {"smallest", "word"}:
            raise ValueError("witness must be 'smallest' or 'word'")
        words = self.orbit_words(obj, ordered=ordered)
        image = min(words)
        if witness == "smallest":
            element = self.transporter(obj, image, ordered=ordered)
        else:
            element = self.materialize(words[image])
        permutation = None
        if isinstance(obj, tuple) and not ordered:
            # where each input position lands in the sorted image
            raw = list(self.apply(element, obj, ordered=True))
            free = list(range(len(image)))
            permutation = tuple(
                free.pop(next(k for k, j in enumerate(free) if image[j] == item))
                for item in raw
            )
        return CanonicalForm(image, element, permutation)

    def orbit_representatives(self, objects: Iterable, *, ordered: bool = False) -> dict:
        """Group ``objects`` by orbit: canonical representative -> members.

        Every orbit is traversed once, so the cost is proportional to the
        total size of the orbits met, not to the group order.
        """

        canonical_of: dict = {}
        buckets: dict = {}
        for obj in objects:
            key = _prepare(obj, ordered)
            representative = canonical_of.get(key)
            if representative is None:
                orbit = self.orbit_words(key, ordered=ordered)
                representative = min(orbit)
                for member in orbit:
                    canonical_of[member] = representative
            buckets.setdefault(representative, []).append(obj)
        return buckets

    # construction

    @classmethod
    def closure(
        cls,
        parent,
        generators: Iterable[Callable[[Element], Element]],
    ) -> AutomorphismGroup[Element, Automorphism[Element]]:
        """The group generated by term maps of a finite parent.

        Elements are image tables; the closure is a breadth-first search on
        the group itself, so the parent's automorphism group must be finite
        (it is, for a finite parent).
        """

        data = _ParentTable(parent)
        generator_tables = list(dict.fromkeys(_table_of(g, data) for g in generators))
        identity = tuple(range(len(data.elements)))
        seen = {identity}
        pending = deque([identity])
        while pending:
            current = pending.popleft()
            for table in generator_tables:
                image = tuple(table[i] for i in current)
                if image not in seen:
                    seen.add(image)
                    pending.append(image)
        elements = [Automorphism(table, data) for table in seen]
        return cls(
            elements,
            apply_term=lambda auto, term: auto(term),
            compose=lambda a, b: a.compose(b),
            inverse=lambda a: a.inverse(),
            identity=Automorphism(identity, data),
            generators=[Automorphism(t, data) for t in generator_tables if t != identity],
            key=(lambda a: (a.generator_images(), a.table)) if data.generator_indices else None,
        )


def _is_identity(element) -> bool:
    return isinstance(element, Automorphism) and element.is_identity()


def _prepare(obj, ordered: bool):
    """Validate an object of the action; sort tuples unless ``ordered``."""

    if isinstance(obj, AdditiveSequence):
        return obj
    if isinstance(obj, tuple) and all(isinstance(item, AdditiveSequence) for item in obj):
        if any(item.parent() is not obj[0].parent() for item in obj):
            raise TypeError("all sequences of a tuple must share one space")
        return obj if ordered else tuple(sorted(obj))
    raise TypeError("expected an additive sequence or a tuple of additive sequences")


# Resolution


def automorphism_group(space: AdditiveSequenceSpace[Element]) -> AutomorphismGroup:
    """Resolve and cache the materialized automorphism group of a space.

    Resolution order: a group supplied by the additive parent
    (``parent.automorphism_group()``), then the closure of the generators of
    the parent's automorphism action (which requires a finite iterable
    parent), otherwise :class:`AutomorphismGroupUnavailable`.
    """

    cached = getattr(space, "_automorphism_group", None)
    if cached is not None:
        return cached
    parent = space.base_parent
    group = None
    provider = getattr(parent, "automorphism_group", None)
    if callable(provider):
        try:
            group = provider()
        except (AttributeError, NotImplementedError):
            group = None
        if group is not None and not isinstance(group, AutomorphismGroup):
            raise TypeError("automorphism_group() must return an AutomorphismGroup")
    if group is None:
        from .orbits import AutomorphismActionUnavailable, automorphism_action

        try:
            action = automorphism_action(space)
        except AutomorphismActionUnavailable:
            action = None
        is_finite = getattr(parent, "is_finite", None)
        if action is not None and callable(is_finite) and is_finite():
            group = AutomorphismGroup.closure(
                parent, [partial(action.apply_term, i) for i in range(len(action))]
            )
    if group is None:
        raise AutomorphismGroupUnavailable(
            "could not materialize an automorphism group for the additive parent"
        )
    space._automorphism_group = group
    return group
