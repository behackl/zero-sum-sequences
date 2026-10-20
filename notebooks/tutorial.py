# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "marimo",
#     "matplotlib",
#     "networkx",
#     "zero-sum-sequences",
# ]
# ///

import marimo

__generated_with = "0.25.2"
app = marimo.App(width="medium", app_title="zero-sum-sequences tutorial")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Tutorial: additive and zero-sum sequences

    `zero_sum_sequences` represents finite multisets over additive groups and
    computes zero-sum atoms, factorizations and their lengths, automorphism
    orbits, and zero-sum constants. The examples use $C_3$,
    $C_2\oplus C_4$, and a four-point support in $C_3^3$.

    This is a [marimo](https://marimo.io) notebook: cells rerun automatically
    when the definitions they use change. It runs in ordinary Python, or in
    the browser from the documentation.
    """)
    return


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    def facts(*rows):
        """Show ``(label, value)`` pairs as a two-column table."""

        lines = ["| | |", "|:--|:--|"]
        lines += [f"| {label} | `{value}` |" for label, value in rows]
        return mo.md("\n".join(lines))

    def format_factorization(factors):
        return "  *  ".join(str(factor) for factor in factors) or "1"

    return facts, format_factorization


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 1. Configure a sequence space

    A finite abelian group can be given by its elements and its addition.
    An `AdditiveSequenceSpace` over it constructs immutable finite multisets
    of group elements.

    Enumerating atoms needs an upper bound for the Davenport constant
    $D(G)$, the maximal length of an atom. When none is supplied, the space
    derives one without an exhaustive search: here $C_3$ is recognized as
    cyclic, so $D(C_3)=3$. An explicit `davenport_bound=` is used as given;
    a bound below $D(G)$ can make results incomplete.
    """)
    return


@app.cell
def _():
    from zero_sum_sequences import (
        AdditiveSequenceSpace,
        AtomCatalogue,
        AutomorphismAction,
        FactorizationRelation,
        FactorizationSolver,
        FiniteAdditiveGroup,
        davenport_constant,
        erdos_ginzburg_ziv_constant,
        eta_constant,
        zero_sum_constant,
    )

    return (
        AdditiveSequenceSpace,
        AtomCatalogue,
        AutomorphismAction,
        FactorizationRelation,
        FactorizationSolver,
        FiniteAdditiveGroup,
        davenport_constant,
        erdos_ginzburg_ziv_constant,
        eta_constant,
        zero_sum_constant,
    )


@app.cell
def _(AdditiveSequenceSpace, FiniteAdditiveGroup, facts):
    modulo_three = FiniteAdditiveGroup(
        range(3),
        zero=0,
        add=lambda left, right: (left + right) % 3,
        coerce=lambda value: int(value) % 3,
    )
    C3 = AdditiveSequenceSpace(modulo_three)
    facts(("group order", modulo_three.order()), ("Davenport bound", C3.davenport_bound))
    return (C3,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 2. Construct and inspect sequences

    Terms are coerced through the group and stored as a canonical
    multiplicity table, so their input order is irrelevant. A nonempty
    zero-sum sequence is an *atom* if it has no proper nonempty zero-sum
    subsequence. `map_terms(mapping)` applies a map term by term.
    """)
    return


@app.cell
def _(C3, facts):
    sample = C3([2, 1, 2, 1])
    negated_sample = sample.map_terms(lambda term: -term)

    facts(
        ("sequence", sample),
        ("terms", tuple(sample)),
        ("length", len(sample)),
        ("support", sample.support),
        ("multiplicities", sample.multiplicities),
        ("sum", sample.total()),
        ("zero-sum", sample.is_zero_sum()),
        ("atom", sample.is_atom()),
        ("termwise negation", negated_sample),
        ("negation fixes the sequence", negated_sample == sample),
    )
    return (sample,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 3. Multiset arithmetic

    Addition combines multiplicities, subtraction removes a subsequence, and
    multiplication by a nonnegative integer repeats a sequence. `divides`
    tests the subsequence relation. Results stay in the same sequence space.
    """)
    return


@app.cell
def _(C3, facts, sample):
    pair = C3([1, 2])
    triple_pair = 3 * pair

    facts(
        ("pair", pair),
        ("pair + pair", pair + pair),
        ("sample - pair", sample - pair),
        ("3 * pair", triple_pair),
        ("pair divides sample", pair.divides(sample)),
        ("3 * pair divides sample", triple_pair.divides(sample)),
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 4. Factorizations and length sets

    Consider $S=1^3 2^3$ over $C_3$. It has two visible factorizations into
    atoms:

    $$S=(1^3)(2^3)=(1\cdot2)^3.$$

    Their lengths are 2 and 3. `length_set()` computes the complete set of
    lengths, while `factorization_witnesses()` keeps one deterministic
    witness for each length.
    """)
    return


@app.cell
def _(C3, facts, format_factorization):
    sequence = C3([1, 1, 1, 2, 2, 2])
    one_cubed = C3([1, 1, 1])
    two_cubed = C3([2, 2, 2])
    mixed_pair = C3([1, 2])

    facts(
        ("length set", sorted(sequence.length_set())),
        *(
            (f"length {_length}", format_factorization(_factors))
            for _length, _factors in sequence.factorization_witnesses().items()
        ),
    )
    return mixed_pair, one_cubed, sequence, two_cubed


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    `factorizations()` enumerates every unordered factorization exactly once.
    For larger examples, prefer the length set or one witness per length.

    Each of these convenience methods solves the problem afresh. For several
    questions about one sequence, keep a solver (section 6).
    """)
    return


@app.cell
def _(facts, format_factorization, sequence):
    facts(
        *(
            (f"length {len(_factors)}", format_factorization(_factors))
            for _factors in sorted(sequence.factorizations(), key=len)
        )
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 5. Factorization relations

    A `FactorizationRelation` records an oriented equality between two
    unordered factorizations. It cancels their common factors and computes
    the standard factorization distance. Adjoining one common copy of
    $1\cdot2$ to the equality above gives

    $$(1\cdot2)(1^3)(2^3) = (1\cdot2)^4.$$
    """)
    return


@app.cell
def _(
    FactorizationRelation,
    facts,
    format_factorization,
    mixed_pair,
    one_cubed,
    two_cubed,
):
    relation = FactorizationRelation(
        source=(one_cubed, two_cubed, mixed_pair),
        target=(mixed_pair,) * 4,
    )

    facts(
        ("common factor", format_factorization(relation.common_factor)),
        ("reduced source", format_factorization(relation.reduced_source)),
        ("reduced target", format_factorization(relation.reduced_target)),
        ("distance", relation.distance),
        ("reduced product", relation.reduced().product),
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 6. Reuse a factorization solver

    A `FactorizationSolver` keeps its discovered atoms and its graph of
    remainders between queries. Exact-length and minimum-length queries use
    targeted searches unless complete length data is already available.
    `statistics` reports the numbers of atom divisors, remainder states, and
    transitions.
    """)
    return


@app.cell
def _(FactorizationSolver, facts, sequence):
    solver = FactorizationSolver(sequence)
    _stats = solver.statistics

    facts(
        ("has length 3", solver.has_factorization_of_length(3)),
        ("minimum length", solver.minimum_factorization_length()),
        ("maximum length", solver.maximum_factorization_length()),
        ("length set", sorted(solver.length_set())),
        ("witness lengths", sorted(solver.factorization_witnesses())),
        ("candidate atoms", _stats.candidate_atoms),
        ("remainder states", _stats.states),
        ("transitions", _stats.transitions),
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 6.1. Ask for just one witness

    `factorization_witness(k)` finds one factorization with exactly `k`
    atoms, or returns `None`. Unlike `factorization_witnesses()`, it does not
    first compute the complete length set.

    For $T=S^3=1^9 2^9$, ask for length seven and check the witness by
    multiplying its factors. A fresh solver ensures that no complete length
    data has been computed yet.
    """)
    return


@app.cell
def _(C3, FactorizationSolver, facts, format_factorization, sequence):
    target = 3 * sequence
    targeted_solver = FactorizationSolver(target)
    witness = targeted_solver.factorization_witness(7)
    assert witness is not None and len(witness) == 7
    assert sum(witness, C3()) == target
    assert all(_atom.is_atom() for _atom in witness)

    facts(
        ("length-7 witness", format_factorization(witness)),
        ("correct product", sum(witness, C3()) == target),
        ("length-5 witness", targeted_solver.factorization_witness(5)),
    )
    return target, targeted_solver


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    A witness can also be required to contain a minimum number of factors
    satisfying a predicate. Here we request at least four cubes (atoms of
    term length three). This is possible at length seven, but not at length
    eight, although length eight itself occurs.

    So `None` from a constrained query means that **both conditions cannot
    hold together**, not necessarily that the length is absent. A positive
    witness can be checked directly; a negative answer relies on an
    exhaustive search and a sufficient bound or catalogue.

    `statistics` and `digraph()` construct the full remainder graph; avoid
    them if only a targeted answer is needed.
    """)
    return


@app.cell
def _(C3, facts, format_factorization, target, targeted_solver):
    def _is_cube(atom):
        return len(atom) == 3

    with_four_cubes = targeted_solver.factorization_witness(
        7, minimum_matching_factors=4, factor_predicate=_is_cube
    )
    assert with_four_cubes is not None and len(with_four_cubes) == 7
    assert sum(with_four_cubes, C3()) == target
    assert sum(map(_is_cube, with_four_cubes)) >= 4

    impossible_profile = targeted_solver.factorization_witness(
        8, minimum_matching_factors=4, factor_predicate=_is_cube
    )
    assert impossible_profile is None
    assert targeted_solver.has_factorization_of_length(8)

    facts(
        ("length 7 with four cubes", format_factorization(with_four_cubes)),
        ("length 8 exists", targeted_solver.has_factorization_of_length(8)),
        ("length 8 with four cubes", impossible_profile),
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 7. Inspect the factorization graph

    `factorization_digraph()` returns the remainder graph as a NetworkX
    directed acyclic graph. Each vertex is a remainder sequence, and each edge
    records the removed atom. A path from the input to the empty sequence
    encodes a factorization.
    """)
    return


@app.cell
def _(C3, facts, sequence):
    import networkx as nx

    dag = sequence.factorization_digraph()
    path_lengths = {
        len(_path) - 1 for _path in nx.all_simple_paths(dag, sequence, C3())
    }

    facts(
        ("vertices", dag.number_of_nodes()),
        ("edges", dag.number_of_edges()),
        ("acyclic", nx.is_directed_acyclic_graph(dag)),
        ("path lengths", sorted(path_lengths)),
        *(
            (f"{_source} → {_target}", f"remove {_atom}")
            for _source, _target, _atom in sorted(
                (str(_s), str(_t), str(_a)) for _s, _t, _a in dag.edges(data="atom")
            )
        ),
    )
    return dag, nx


@app.cell
def _(dag, nx):
    import matplotlib.pyplot as plt

    _positions = {}
    for _depth, _nodes in enumerate(nx.topological_generations(dag)):
        _ordered = sorted(_nodes, key=str)
        _offset = (len(_ordered) - 1) / 2
        for _index, _node in enumerate(_ordered):
            _positions[_node] = (_index - _offset, -_depth)

    _figure, _axis = plt.subplots(figsize=(7, 4.5))
    nx.draw_networkx_edges(
        dag, pos=_positions, node_size=900, node_shape="s", arrowsize=18, ax=_axis
    )
    nx.draw_networkx_labels(
        dag,
        pos=_positions,
        labels={_node: str(_node) for _node in dag},
        font_size=10,
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "#e8eef7", "edgecolor": "#607089"},
        ax=_axis,
    )
    nx.draw_networkx_edge_labels(
        dag,
        pos=_positions,
        edge_labels={(_s, _t): str(_a) for _s, _t, _a in dag.edges(data="atom")},
        font_size=10,
        rotate=False,
        ax=_axis,
    )
    _axis.set_axis_off()
    _figure.tight_layout()
    _axis
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 8. Enumerate atoms

    Every atom over a finite group $G$ has length at most $D(G)$.
    `enumerate_atom_catalogue()` lists all atoms up to the space's Davenport
    bound. Like the factorization queries, it omits the singleton of the
    identity, which is an atom on its own.

    For $C_2\oplus C_4$ the derived bound is exact,
    $D(C_2\oplus C_4)=1+(2-1)+(4-1)=5$, and the complete enumeration is
    small.
    """)
    return


@app.cell
def _(AdditiveSequenceSpace, FiniteAdditiveGroup, facts):
    from collections import Counter

    group = FiniteAdditiveGroup.cyclic_product(2, 4)
    C2xC4 = AdditiveSequenceSpace(group)
    c2_c4_catalogue = C2xC4.enumerate_atom_catalogue()

    facts(
        ("Davenport bound", C2xC4.davenport_bound),
        ("atoms", len(c2_c4_catalogue)),
        *(
            (f"atoms of length {_length}", _count)
            for _length, _count in sorted(Counter(map(len, c2_c4_catalogue)).items())
        ),
        (
            "all are atoms without the identity",
            all(_a.is_atom() and group.zero() not in _a for _a in c2_c4_catalogue),
        ),
    )
    return C2xC4, Counter, c2_c4_catalogue, group


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 9. Automorphism orbits

    Automorphisms act term by term. `FiniteAdditiveGroup.cyclic_product(...)`
    supplies generators of the full automorphism group.

    The eight atoms of length 5 over $C_2\oplus C_4$ form one orbit. An orbit
    witness is a shortest word in the generators mapping one atom to another;
    it can show the induced images of the standard additive generators.
    """)
    return


@app.cell
def _(AutomorphismAction, C2xC4, c2_c4_catalogue, facts, group):
    maximal_atom = C2xC4([(0, 1)] * 3 + [(1, 0), (1, 1)])
    other_maximal_atom = C2xC4([(0, 1), (1, 0)] + [(1, 1)] * 3)
    maximal_orbit = maximal_atom.orbit()
    orbit_witness = maximal_atom.orbit_witness(other_maximal_atom)
    action = AutomorphismAction(group.automorphism_generators())

    facts(
        ("orbit size", len(maximal_orbit)),
        ("same orbit", maximal_atom.is_in_same_orbit(other_maximal_atom)),
        ("generator word", orbit_witness.generator_indices),
        *(
            (f"generator {_g}", f"↦ {_image}")
            for _g, _image in orbit_witness.additive_generator_images()
        ),
        (
            "orbit is every length-5 atom",
            set(maximal_orbit) == {_a for _a in c2_c4_catalogue if len(_a) == 5},
        ),
        (
            "witness maps atom to target",
            action.apply_word(maximal_atom, orbit_witness) == other_maximal_atom,
        ),
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 10. Reuse a complete atom catalogue

    An `AtomCatalogue` indexes precomputed atoms for reuse. For the earlier
    $C_3$ input, the relevant atoms are $1^3$, $2^3$, and $1\cdot2$.
    Complete factorization results require every atom divisor of the
    sequence; the catalogue's completeness is not checked.
    """)
    return


@app.cell
def _(
    AtomCatalogue,
    C3,
    FactorizationSolver,
    facts,
    mixed_pair,
    one_cubed,
    sequence,
    two_cubed,
):
    catalogue = AtomCatalogue(C3, (one_cubed, two_cubed, mixed_pair))
    catalogue_solver = FactorizationSolver(sequence, atom_catalogue=catalogue)

    facts(
        ("catalogue atoms", [str(_a) for _a in catalogue]),
        ("divisors of S", [str(_a) for _a in catalogue.divisors(sequence)]),
        ("length set", sorted(sequence.length_set(atom_catalogue=catalogue))),
        ("solver statistics", catalogue_solver.statistics),
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 11. Atoms on a prescribed support

    Pass `support=G0` to enumerate only atoms whose terms belong to a subset
    $G_0$ of the group, without enumerating and filtering the full
    catalogue. **The support need not be a subgroup:** sums and inverses are
    still taken in the ambient group.

    Take $G=C_3^3$, with $D(G)=7$, and
    $G_0=\{e_1,e_2,e_3,e_1+e_2+e_3\}$. There are just six atoms on this
    support: four cubes, one atom of term length five, and one of term
    length seven.
    """)
    return


@app.cell
def _(AdditiveSequenceSpace, Counter, FiniteAdditiveGroup, facts):
    C3cubed = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(3, 3, 3))
    e1, e2, e3, e4 = (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)
    G0 = [e1, e2, e3, e4]
    support_catalogue = C3cubed.enumerate_atom_catalogue(support=G0)

    assert Counter(map(len, support_catalogue)) == {3: 4, 5: 1, 7: 1}
    assert all(set(_a.support) <= set(G0) for _a in support_catalogue)
    facts(*((f"term length {len(_a)}", _a) for _a in support_catalogue))
    return C3cubed, e1, e2, e3, e4, support_catalogue


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    The restricted catalogue is complete for **every sequence on this
    support**, given a sufficient bound: every factor of such a sequence
    also has its terms in $G_0$. It is generally not complete for sequences
    with terms outside $G_0$.

    For the length-seven atom $U$, the sequence $U^3$ has length set
    $\{3,5,7\}$. We first ask for a length-five witness and for the missing
    length four, and compute the complete length set only afterwards.
    """)
    return


@app.cell
def _(
    C3cubed,
    FactorizationSolver,
    e1,
    e2,
    e3,
    e4,
    facts,
    format_factorization,
    support_catalogue,
):
    U = C3cubed([e1, e1, e2, e2, e3, e3, e4])
    support_solver = FactorizationSolver(3 * U, atom_catalogue=support_catalogue)
    middle_witness = support_solver.factorization_witness(5)
    assert middle_witness is not None and len(middle_witness) == 5
    assert sum(middle_witness, C3cubed()) == 3 * U
    assert support_solver.factorization_witness(4) is None
    assert support_solver.length_set() == {3, 5, 7}

    facts(
        ("length-5 witness", format_factorization(middle_witness)),
        ("length-4 witness", support_solver.factorization_witness(4)),
        ("complete length set", sorted(support_solver.length_set())),
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 12. Zero-sum constants

    The package computes classical zero-sum constants exactly:

    - $D(G)$, the least length forcing a nonempty zero-sum subsequence;
    - $\eta(G)$, the least length forcing one of length at most $\exp(G)$;
    - $\mathsf s(G)$, the least length forcing one of length exactly
      $\exp(G)$ (the Erdős–Ginzburg–Ziv constant).

    Known formulas are used where they apply; otherwise the functions search
    exhaustively, reduced by automorphisms. Exhaustive searches grow quickly
    with the group. For a rank-two group $C_m\oplus C_n$ with $m\mid n$,
    the known values are $D=m+n-1$, $\eta=2m+n-2$, and $\mathsf s=2m+2n-3$.

    Choose a group:
    """)
    return


@app.cell
def _(mo):
    group_choice = mo.ui.dropdown(
        options={
            "C₂ ⊕ C₄": (2, 4),
            "C₃ ⊕ C₃": (3, 3),
            "C₂ ⊕ C₆": (2, 6),
            "C₄ ⊕ C₄": (4, 4),
            "C₂ ⊕ C₂ ⊕ C₂": (2, 2, 2),
        },
        value="C₂ ⊕ C₄",
        label="group",
    )
    group_choice
    return (group_choice,)


@app.cell
def _(
    AdditiveSequenceSpace,
    FiniteAdditiveGroup,
    davenport_constant,
    erdos_ginzburg_ziv_constant,
    eta_constant,
    facts,
    group_choice,
):
    _moduli = group_choice.value
    chosen_space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(*_moduli))
    _rows = [
        ("D(G)", davenport_constant(chosen_space)),
        ("η(G)", eta_constant(chosen_space)),
        ("s(G)", erdos_ginzburg_ziv_constant(chosen_space)),
    ]
    if len(_moduli) == 2:
        _m, _n = _moduli
        _rows.append(("rank-two formulas", (_m + _n - 1, 2 * _m + _n - 2, 2 * _m + 2 * _n - 3)))
    facts(*_rows)
    return (chosen_space,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    `zero_sum_constant(space, L)` is the general form: the least length
    forcing a nonempty zero-sum subsequence whose length lies in $L$. Its
    result also lists the extremal sequences, one of length one less than the
    constant per automorphism orbit. Here $L=\{\exp(G)\}$ gives
    $\mathsf s(G)$ again:
    """)
    return


@app.cell
def _(chosen_space, facts, zero_sum_constant):
    _exponent = chosen_space.base_parent.exponent()
    egz = zero_sum_constant(chosen_space, {_exponent})

    facts(
        ("allowed lengths", sorted(egz.lengths)),
        ("constant", egz.value),
        ("extremal orbits", len(egz.extremal)),
        ("first extremal sequence", egz.extremal[0]),
    )
    return


if __name__ == "__main__":
    app.run()
