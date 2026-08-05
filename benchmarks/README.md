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
sage -python -m benchmarks.benchmark_factorization --enumerate
```
