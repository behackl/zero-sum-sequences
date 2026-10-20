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
    orbits, and zero-sum constants. The examples use $C_3$, $C_2\oplus C_4$,
    a four-point support in $C_3^3$, and finally a few small groups for
    zero-sum constants.

    This is a [marimo](https://marimo.io) notebook: when you edit a cell, the
    cells that depend on it rerun automatically. It runs in ordinary Python,
    or in the browser from the documentation.
    """)
    return


@app.cell
def _():
    import math

    import marimo as mo

    from zero_sum_sequences import (
        AdditiveSequenceSpace,
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
        AutomorphismAction,
        FactorizationRelation,
        FactorizationSolver,
        FiniteAdditiveGroup,
        davenport_constant,
        erdos_ginzburg_ziv_constant,
        eta_constant,
        math,
        mo,
        zero_sum_constant,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 1. Configure a sequence space

    A finite abelian group can be given by its elements and its addition.
    An `AdditiveSequenceSpace` over it constructs immutable finite multisets
    of group elements.

    Enumerating atoms needs an upper bound for the Davenport constant
    $D(G)$, the maximal length of an atom. When none is supplied, the space
    derives one without an exhaustive search; here $C_3$ is recognized as
    cyclic, so the bound is $D(C_3)=3$. An explicit `davenport_bound=` is used
    as given, and a bound below $D(G)$ can make results incomplete.
    """)
    return


@app.cell
def _(AdditiveSequenceSpace, FiniteAdditiveGroup):
    modulo_three = FiniteAdditiveGroup(
        range(3),
        zero=0,
        add=lambda left, right: (left + right) % 3,
        coerce=lambda value: int(value) % 3,
    )
    C3 = AdditiveSequenceSpace(modulo_three)
    print("group order:", modulo_three.order())
    print("Davenport bound:", C3.davenport_bound)
    return (C3,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 2. Construct and inspect sequences

    Terms are coerced through the group and stored as a canonical
    multiplicity table, so their input order is irrelevant. A nonempty
    zero-sum sequence is an *atom* if it has no proper nonempty zero-sum
    subsequence. `-S` negates every term; more generally,
    `map_terms(mapping)` applies a map term by term.
    """)
    return


@app.cell
def _(C3):
    sample = C3(2, 1, 2, 1)
    negated_sample = -sample

    print("sequence:", sample)
    print("parent is C3:", sample.parent() is C3)
    print("terms:", tuple(sample))
    print("length:", len(sample))
    print("support:", sample.support)
    print("multiplicities:", sample.multiplicities)
    print("sum:", sample.sum())
    print("zero-sum:", sample.is_zero_sum())
    print("atom:", sample.is_atom())
    print("termwise negation:", negated_sample)
    print("negation fixes the sequence:", negated_sample == sample)
    return (sample,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 3. Products of sequences

    Sequences are written multiplicatively: `S * T` is the product of two
    sequences, `S ** k` the $k$-th power, and `S / T` is $S T^{-1}$ for a
    divisor $T$ of $S$. `T.divides(S)` tests whether $T$ is a subsequence of
    $S$. Results stay in the same sequence space.
    """)
    return


@app.cell
def _(C3, sample):
    pair = C3(1, 2)
    pair_cubed = pair ** 3

    print("pair:", pair)
    print("pair * pair:", pair * pair)
    print("sample / pair:", sample / pair)
    print("pair ** 3:", pair_cubed)
    print("pair divides sample:", pair.divides(sample))
    print("pair ** 3 divides sample:", pair_cubed.divides(sample))
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
def _(C3):
    one_cubed = C3(1) ** 3
    two_cubed = C3(2) ** 3
    mixed_pair = C3(1, 2)
    sequence = one_cubed * two_cubed


    def format_factorization(factors):
        return "  *  ".join(str(factor) for factor in factors) or "1"


    print("length set:", sorted(sequence.length_set()))
    for _length, _factors in sequence.factorization_witnesses().items():
        print(f"length {_length}: {format_factorization(_factors)}")
    return format_factorization, mixed_pair, one_cubed, sequence, two_cubed


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    `factorizations()` enumerates every unordered factorization exactly once.
    Here these are just the two factorizations above. For larger examples,
    prefer the length set or one witness per length.
    """)
    return


@app.cell
def _(format_factorization, sequence):
    for _factors in sorted(sequence.factorizations(), key=len):
        print(f"length {len(_factors)}: {format_factorization(_factors)}")
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
    format_factorization,
    mixed_pair,
    one_cubed,
    two_cubed,
):
    relation = FactorizationRelation(
        source=(one_cubed, two_cubed, mixed_pair),
        target=(mixed_pair,) * 4,
    )

    print("common factor:", format_factorization(relation.common_factor))
    print("reduced source:", format_factorization(relation.reduced_source))
    print("reduced target:", format_factorization(relation.reduced_target))
    print("distance:", relation.distance)
    print("reduced product:", relation.reduced().product)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 6. Reuse a factorization solver

    Each of the sequence methods above solves its problem afresh. For several
    questions about one sequence, keep a `FactorizationSolver`: it retains
    its discovered atoms and its graph of remainders between queries.
    Exact-length and minimum-length queries use targeted searches unless
    complete length data is already available. `statistics` reports the
    numbers of atom divisors, remainder states, and transitions.
    """)
    return


@app.cell
def _(FactorizationSolver, sequence):
    solver = FactorizationSolver(sequence)
    stats = solver.statistics

    print("has length 3:", solver.has_factorization_of_length(3))
    print("minimum length:", solver.minimum_factorization_length())
    print("maximum length:", solver.maximum_factorization_length())
    print("length set:", sorted(solver.length_set()))
    print("witness lengths:", sorted(solver.factorization_witnesses()))
    print("candidate atoms:", stats.candidate_atoms)
    print("remainder states:", stats.states)
    print("transitions:", stats.transitions)
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
    data has been computed yet. Since $T$ only factors into 6 to 9 atoms, a
    request for length five returns `None`.
    """)
    return


@app.cell
def _(FactorizationSolver, format_factorization, math, sequence):
    target = sequence ** 3
    targeted_solver = FactorizationSolver(target)
    witness = targeted_solver.factorization_witness(7)
    assert witness is not None and len(witness) == 7
    assert math.prod(witness) == target
    assert all(atom.is_atom() for atom in witness)

    print("length-7 witness:", format_factorization(witness))
    print("correct product:", math.prod(witness) == target)
    print("length-5 witness:", targeted_solver.factorization_witness(5))
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
def _(format_factorization, math, target, targeted_solver):
    def is_cube(atom):
        return len(atom) == 3


    with_four_cubes = targeted_solver.factorization_witness(
        7, minimum_matching_factors=4, factor_predicate=is_cube
    )
    assert with_four_cubes is not None and len(with_four_cubes) == 7
    assert math.prod(with_four_cubes) == target
    assert sum(map(is_cube, with_four_cubes)) >= 4

    impossible_profile = targeted_solver.factorization_witness(
        8, minimum_matching_factors=4, factor_predicate=is_cube
    )
    assert impossible_profile is None
    assert targeted_solver.has_factorization_of_length(8)

    print("length 7 with four cubes:", format_factorization(with_four_cubes))
    print("length 8 exists:", targeted_solver.has_factorization_of_length(8))
    print("length 8 with four cubes:", impossible_profile)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 7. Inspect the factorization graph

    `factorization_digraph()` returns the remainder graph as a NetworkX
    directed acyclic graph. Each vertex is a remainder sequence, and each
    edge records the removed atom. A path from the input to the empty
    sequence encodes a factorization.
    """)
    return


@app.cell
def _(C3, sequence):
    import networkx as nx

    dag = sequence.factorization_digraph()
    edges = sorted(
        (str(source), str(end), str(atom))
        for source, end, atom in dag.edges(data="atom")
    )

    print(f"vertices: {dag.number_of_nodes()}; edges: {dag.number_of_edges()}")
    for _source, _end, _atom in edges:
        print(f"{_source}  -- remove {_atom} -->  {_end}")

    path_lengths = {len(path) - 1 for path in nx.all_simple_paths(dag, sequence, C3())}
    print("acyclic:", nx.is_directed_acyclic_graph(dag))
    print("path lengths:", sorted(path_lengths))
    return dag, nx


@app.cell
def _(dag, nx):
    import matplotlib.pyplot as plt

    positions = {}
    for depth, nodes in enumerate(nx.topological_generations(dag)):
        ordered_nodes = sorted(nodes, key=str)
        offset = (len(ordered_nodes) - 1) / 2
        for index, node in enumerate(ordered_nodes):
            positions[node] = (index - offset, -depth)

    figure, axis = plt.subplots(figsize=(7, 4.5))
    nx.draw_networkx_edges(
        dag, pos=positions, node_size=900, node_shape="s", arrowsize=18, ax=axis
    )
    nx.draw_networkx_labels(
        dag,
        pos=positions,
        labels={node: str(node) for node in dag},
        font_size=10,
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "#e8eef7", "edgecolor": "#607089"},
        ax=axis,
    )
    nx.draw_networkx_edge_labels(
        dag,
        pos=positions,
        edge_labels={(source, end): str(atom) for source, end, atom in dag.edges(data="atom")},
        font_size=10,
        rotate=False,
        ax=axis,
    )
    axis.set_axis_off()
    figure.tight_layout()
    axis
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
def _(AdditiveSequenceSpace, FiniteAdditiveGroup):
    from collections import Counter

    group = FiniteAdditiveGroup.cyclic_product(2, 4)
    C2xC4 = AdditiveSequenceSpace(group)
    c2_c4_catalogue = C2xC4.enumerate_atom_catalogue()

    print("Davenport bound:", C2xC4.davenport_bound)
    print("atoms:", len(c2_c4_catalogue))
    for _length, _count in sorted(Counter(map(len, c2_c4_catalogue)).items()):
        print(f"length {_length}: {_count}")
    print(
        "all are atoms without the identity:",
        all(atom.is_atom() and group.zero() not in atom for atom in c2_c4_catalogue),
    )
    return C2xC4, Counter, c2_c4_catalogue, group


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 9. Reuse a complete atom catalogue

    An `AtomCatalogue` indexes atoms for reuse: a solver given a catalogue
    looks up atom divisors instead of searching for them. Results are
    complete when the catalogue contains every atom divisor of the input,
    as the full catalogue above does; this is not checked.

    Take the length-five atom $A=(0,1)^3(1,0)(1,1)$ and its negative
    $B=-A$. Their product has factorizations of lengths 2, 4, and 5, but not 3:
    a set of lengths need not be an interval.
    """)
    return


@app.cell
def _(C2xC4, FactorizationSolver, c2_c4_catalogue, format_factorization):
    A = C2xC4((0, 1)) ** 3 * C2xC4((1, 0), (1, 1))
    B = -A
    product = A * B
    catalogue_solver = FactorizationSolver(product, atom_catalogue=c2_c4_catalogue)

    print("A:", A)
    print("B:", B)
    print("atom divisors in the catalogue:", len(list(c2_c4_catalogue.divisors(product))))
    print("length set:", sorted(catalogue_solver.length_set()))
    for _length, _factors in catalogue_solver.factorization_witnesses().items():
        print(f"length {_length}: {format_factorization(_factors)}")
    return (A,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## 10. Automorphism orbits

    Automorphisms act term by term. `FiniteAdditiveGroup.cyclic_product(...)`
    supplies generators of the full automorphism group.

    The eight atoms of length 5 over $C_2\oplus C_4$, including $A$, form one
    orbit. An orbit witness is a shortest word in the generators mapping one
    atom to another; `show()` prints the induced images of the standard
    additive generators.
    """)
    return


@app.cell
def _(A, AutomorphismAction, C2xC4, c2_c4_catalogue, group):
    other_maximal_atom = C2xC4((0, 1), (1, 0)) * C2xC4((1, 1)) ** 3
    maximal_orbit = A.orbit()
    orbit_witness = A.orbit_witness(other_maximal_atom)

    print("orbit size:", len(maximal_orbit))
    print("same orbit:", A.is_in_same_orbit(other_maximal_atom))
    print("generator word:", orbit_witness.generator_indices)
    orbit_witness.show()

    action = AutomorphismAction(group.automorphism_generators())
    length_five_atoms = {atom for atom in c2_c4_catalogue if len(atom) == 5}
    print("orbit contains every length-5 atom:", set(maximal_orbit) == length_five_atoms)
    print(
        "witness maps A to the target:",
        action.apply_word(A, orbit_witness) == other_maximal_atom,
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
def _(AdditiveSequenceSpace, Counter, FiniteAdditiveGroup):
    C3cubed = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(3, 3, 3))
    e1, e2, e3, e4 = (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)
    G0 = [e1, e2, e3, e4]
    support_catalogue = C3cubed.enumerate_atom_catalogue(support=G0)

    assert Counter(map(len, support_catalogue)) == {3: 4, 5: 1, 7: 1}
    assert all(set(atom.support) <= set(G0) for atom in support_catalogue)
    for atom in support_catalogue:
        print(f"term length {len(atom)}: {atom}")
    return C3cubed, e1, e2, e3, e4, support_catalogue


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    The restricted catalogue is complete for **every sequence on this
    support**, given a sufficient bound: every factor of such a sequence also
    has its terms in $G_0$. It is generally not complete for sequences with
    terms outside $G_0$.

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
    format_factorization,
    math,
    support_catalogue,
):
    U = C3cubed(e1, e1, e2, e2, e3, e3, e4)
    support_solver = FactorizationSolver(U ** 3, atom_catalogue=support_catalogue)
    middle_witness = support_solver.factorization_witness(5)
    assert middle_witness is not None and len(middle_witness) == 5
    assert math.prod(middle_witness) == U ** 3
    assert support_solver.factorization_witness(4) is None
    assert support_solver.length_set() == {3, 5, 7}

    print("length-5 witness:", format_factorization(middle_witness))
    print("length-4 witness:", support_solver.factorization_witness(4))
    print("complete length set:", sorted(support_solver.length_set()))
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
    exhaustively, reduced by automorphisms. To try another group, edit
    `moduli` below. Keep the group small: exhaustive searches grow quickly
    with its order.
    """)
    return


@app.cell
def _(
    AdditiveSequenceSpace,
    FiniteAdditiveGroup,
    davenport_constant,
    erdos_ginzburg_ziv_constant,
    eta_constant,
):
    moduli = (2, 4)
    chosen_space = AdditiveSequenceSpace(FiniteAdditiveGroup.cyclic_product(*moduli))

    print("D(G):", davenport_constant(chosen_space))
    print("eta(G):", eta_constant(chosen_space))
    print("s(G):", erdos_ginzburg_ziv_constant(chosen_space))
    return (chosen_space,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    `zero_sum_constant(space, L)` is the general form: the least length
    forcing a nonempty zero-sum subsequence whose length lies in $L$. Its
    result also lists the extremal sequences, of length one less than the
    constant, one per automorphism orbit. With $L=\{\exp(G)\}$ it gives
    $\mathsf s(G)$ again. Unlike factorizations, these searches allow the
    identity as a term, and extremal sequences can contain it.
    """)
    return


@app.cell
def _(chosen_space, zero_sum_constant):
    egz = zero_sum_constant(chosen_space, {chosen_space.base_parent.exponent()})

    print("allowed lengths:", sorted(egz.lengths))
    print("constant:", egz.value)
    print("extremal orbits:", len(egz.extremal))
    print("first extremal sequence:", egz.extremal[0])
    return


if __name__ == "__main__":
    app.run()
