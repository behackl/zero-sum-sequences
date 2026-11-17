Atoms and factorizations
========================

.. currentmodule:: zero_sum_sequences

Atom catalogues, factorization queries and relations between
factorizations. The factorization methods of :class:`AdditiveSequence`
build a new :class:`FactorizationSolver` on every call.

.. autosummary::
   :nosignatures:

   FactorizationSolver
   FactorizationCache
   AtomCatalogue
   FactorizationRelation

.. autoclass:: FactorizationSolver
   :members:

.. autoclass:: FactorizationStatistics
   :members:

.. autoclass:: FactorizationCache
   :members:

.. autoclass:: FactorizationCacheStatistics
   :members:

.. autoclass:: AnchoredCanonicalizer
   :members:

.. autoclass:: AtomCatalogue
   :members:

.. autoclass:: AtomOrbit
   :members:

.. autoclass:: FactorizationRelation
   :members:

.. py:data:: Factorization

   A factorization: a tuple of atoms whose product is the factorized
   sequence. Returned by the witness and enumeration queries.
