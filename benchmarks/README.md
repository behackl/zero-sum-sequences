# Factorization benchmark corpus

The benchmark corpus is shared by the test suite and timing runner. It
separates input size from combinatorial difficulty and checks the sole
supported memoized factorization implementation.

| Case | Tier | Terms | Expected length set |
|---|---|---:|---|
| `c3-atom` | short | 3 | `{1}` |
| `c3-balanced-block` | short | 6 | `{2,3}` |
| `rank-three-inverse-pair` | short | 14 | `{2,3,4,5,7}` |
| `c3-balanced-power-7` | long | 42 | `[14,21]` |
| `c3-balanced-power-40` | long | 240 | `[80,120]` |
| `rank-three-inverse-power-3` | long | 42 | `[6,21]` |
| `c3-balanced-power-250` | very long | 1,500 | `[500,750]` |
| `c3-pure-power-5000` | very long | 15,000 | `{5000}` |

For the balanced `C_3` family, write `A = 1^3`, `B = 2^3`, and
`C = 1 * 2`. A factorization of `(A * B)^n` has `x` copies each of `A`
and `B`, and `3(n-x)` copies of `C`, for `0 <= x <= n`. Its length set is
therefore exactly `[2n,3n]`, and it has `n+1` unordered factorizations. The
pure-power case has one factorization of length 5,000.

The rank-three cases exercise a less artificial support and a larger
candidate-atom set. The runner reports candidate indexing, state-graph solve,
and safe exhaustive-enumeration timings while validating all expected results.

Run:

```console
python -m benchmarks.benchmark_factorization --enumerate
```

## Operations benchmark

`benchmark_operations.py` times catalogue enumeration, atom tests, length
sets with and without a catalogue, automorphism operations, and zero-sum
constants. The cases cover cyclic groups, groups of rank two, elementary
abelian groups and `C_3^3` (27 911 atoms, 11 232 automorphisms), with both
generic and highly symmetric inputs.

```console
python -m benchmarks.benchmark_operations          # quick cases
python -m benchmarks.benchmark_operations --large  # adds larger catalogues
```

Reference timings in milliseconds for 0.3.0 and 0.3.1 (one core, Python
3.13, October 2026), from comparable workloads run with both versions:

| Workload | 0.3.0 | 0.3.1 |
|---|---:|---:|
| catalogue `C_12` | 18 000 | 1.9 |
| catalogue `C_2 + C_8` | 11 700 | 9.2 |
| catalogue `C_4 + C_4` | 1 210 | 7.4 |
| catalogue `C_5 + C_5` | > 60 000 | 340 |
| catalogue `C_3^3` | about 28 000 | 320 |
| `is_atom`, 3000 random sequences over `C_5 + C_5` | 33 | 5.4 |
| `is_atom`, high multiplicities over `C_200` | 17 | 11 |
| length sets without catalogue over `C_5 + C_5` | 355 | 74 |
| length set of `1^18 5^18 7^18` over `C_18` | 31 | 6.4 |
| length sets with catalogue over `C_3^3`, 100 products | 3 880 | 725 |
| canonical forms, generic sequences over `C_2^4` | 2 060 | 80 |
| canonical forms, symmetric sequences over `C_2^4` | 1.0 | 0.2 |
| canonical forms, generic sequences over `C_5 + C_5` | 95 | 2.1 |
| canonical forms, atom pairs over `C_3^3` | 2 160 | 29 |
| stabilizers, generic sequences over `C_2^4` | 359 | 64 |
| orbit classes of 300 sequences over `C_5 + C_5` | 958 | 112 |
| `eta(C_3^3)` | 3 780 | 850 |
