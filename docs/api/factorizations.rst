Atoms and factorizations
========================

.. currentmodule:: zero_sum_sequences

The sequence methods such as :meth:`AdditiveSequence.length_set` are
enough for single queries. For several queries about one sequence, keep a
:class:`FactorizationSolver`. For many sequences over one group, enumerate
an :class:`AtomCatalogue` once and answer queries through a
:class:`FactorizationCache`.

.. autosummary::
   :nosignatures:

   FactorizationSolver
   AtomCatalogue
   FactorizationCache
   FactorizationRelation

.. autoclass:: FactorizationSolver
   :members:

.. autoclass:: FactorizationStatistics
   :members:

.. autoclass:: AtomCatalogue
   :members:

.. autoclass:: AtomOrbit
   :members:

.. autoclass:: FactorizationCache
   :members:

.. autoclass:: FactorizationCacheStatistics
   :members:

.. autoclass:: AnchoredCanonicalizer
   :members:

.. autoclass:: FactorizationRelation
   :members:

.. py:data:: Factorization

   A factorization: a tuple of atoms, sorted, whose product is the
   factorized sequence. Returned by the witness and enumeration queries.
