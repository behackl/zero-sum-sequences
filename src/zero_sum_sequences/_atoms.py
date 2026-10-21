"""Atom tests and atom enumeration through zero-sum free sequences.

Both rest on one observation: a zero-sum sequence ``S != 1`` is an atom if and
only if ``S g^{-1}`` is zero-sum free for one (equivalently, every) term ``g``
of ``S``. Indeed, a proper nonempty zero-sum subsequence ``T`` and its
complement ``S T^{-1}`` are both zero-sum, and one of them misses the removed
copy of ``g``.

A sequence is zero-sum free when no nonempty subsequence sums to zero. Its
set of nonempty subsequence sums ``Σ`` extends term by term,
``Σ(T g) = Σ(T) ∪ {g} ∪ (Σ(T) + g)``, and ``T g`` stays zero-sum free exactly
when ``-g ∉ Σ(T)`` and ``g != 0``.
"""

from __future__ import annotations

# Index tables for enumeration are built when |G| * |support| stays below this.
_TABLE_LIMIT = 1 << 20


def is_atom_by_sums(sequence) -> bool:
    """Return whether a nonempty zero-sum ``sequence`` is an atom.

    Needs only the parent's addition and zero: since ``Σ(T)`` never contains
    zero along the way, ``T g`` stays zero-sum free exactly when zero is not
    in the updated set of sums.
    """

    space = sequence.parent()
    add, zero = space._add_canonical(), space._canonical_zero()
    items = list(sequence._items)
    term, count = items[-1]
    items[-1] = (term, count - 1)  # remove one copy of the largest term
    sums: set = set()
    for term, count in items:
        for _ in range(count):
            # Σ(T g) = Σ(T) ∪ {g} ∪ (Σ(T) + g), shifting the old Σ(T).
            sums |= {term} | {add(total, term) for total in sums}
            if zero in sums:
                return False
    return True


def enumerate_atoms(space, terms, allowed, bound) -> list[tuple]:
    """Return the item tuples of all atoms of length at most ``bound`` on ``allowed``.

    ``terms`` are all canonical elements of the finite parent, sorted;
    ``allowed`` is a sorted tuple of nonzero canonical terms. Each atom is
    produced once: as its sorted terms without one copy of the largest,
    which form a zero-sum free prefix, followed by that largest term.
    """

    arithmetic = space._group_arithmetic()
    add, negate = arithmetic.add_canonical, arithmetic.negate
    if not allowed or bound < 2:
        return []
    if len(terms) * len(allowed) <= _TABLE_LIMIT:
        index = {term: position for position, term in enumerate(terms)}
        labels = tuple(index[term] for term in allowed)
        negative = tuple(index[negate(term)] for term in terms)
        columns = {
            label: tuple(index[add(term, terms[label])] for term in terms)
            for label in labels
        }

        def shift(label):
            column = columns[label]
            return lambda sums: {column[total] for total in sums}

        def add_label(total, label):
            return columns[label][total]

        def negative_of(label):
            return negative[label]

        def term_of(label):
            return terms[label]

    else:
        labels = allowed

        def shift(label):
            return lambda sums: {add(total, label) for total in sums}

        add_label = add
        negative_of = negate

        def term_of(label):
            return label

    allowed_labels = frozenset(labels)
    shifts = {label: shift(label) for label in labels}
    found = []
    # Depth-first over sorted zero-sum free prefixes: (prefix, Σ(prefix), σ(prefix)).
    stack = [((label,), {label}, label) for label in reversed(labels)]
    while stack:
        prefix, sums, total = stack.pop()
        last = prefix[-1]
        final = negative_of(total)
        if final in allowed_labels and final >= last:
            found.append((*prefix, final))
        if len(prefix) + 1 >= bound:
            continue
        for label in labels:
            if label < last or negative_of(label) in sums:
                continue
            extended = sums | {label} | shifts[label](sums)
            stack.append(((*prefix, label), extended, add_label(total, label)))

    atoms = []
    for labels_of_atom in found:
        counts: dict = {}
        for label in labels_of_atom:
            counts[label] = counts.get(label, 0) + 1
        atoms.append(tuple((term_of(label), count) for label, count in sorted(counts.items())))
    return atoms
