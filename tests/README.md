# Test suite

Run the ordinary-Python suite from the repository root:

```console
uv run python -m pytest -q
```

With SageMath available, include the native-parent tests:

```console
PYTHONPATH="$PWD/src:$PWD" sage -python -m pytest -q -p no:cacheprovider
```

## Coverage ownership

- **Component suites** own public behavior, serialization, and named regressions
  for sequences, catalogues, arithmetic, invariants, symmetry, and factorization.
  Optional-bound integration remains in `test_optional_bounds.py`.
- **`test_parent_protocols.py`** covers foreign-parent arithmetic and the
  distinction between unavailable capabilities and broken providers.
- **`test_validation.py`** covers shared exact-integer rules. Component tests
  check delegation and API-specific interactions, notably validation before
  cache lookup and rejection of booleans/floats that compare equal to cached keys.
- **`test_coordinate_oracle.py`** independently computes sums and enumerates
  sequences in cyclic coordinates. Its mathematical ground truth does not use
  package sequence predicates. Larger invariant cases complement that oracle.
- **`test_generic_parents.py` and `test_orbits_sage.py`** cover native Sage parents
  and are skipped when Sage is unavailable.

Use representative inputs for shared rules; expand a parameter grid when the
interaction itself matters. Work-limit tests use operation counts or forbidden
operations, not timing thresholds. Private state is inspected only where it
expresses a specific caching or allocation guarantee.

Timing comparisons live in `benchmarks/`; its small correctness
corpus is also exercised by `test_factorization_benchmarks.py`.
