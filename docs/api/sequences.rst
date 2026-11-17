Groups and sequences
====================

.. currentmodule:: zero_sum_sequences

A :class:`FiniteAdditiveGroup` describes the base group; compatible Sage
parents can be used instead. An :class:`AdditiveSequenceSpace` constructs
sequences over a group, and :class:`AdditiveSequence` carries the
per-sequence queries: arithmetic, atom tests, factorizations and orbits.

.. autosummary::
   :nosignatures:

   FiniteAdditiveGroup
   AdditiveSequenceSpace
   AdditiveSequence

.. autoclass:: FiniteAdditiveGroup
   :members:

.. autoclass:: AdditiveSequenceSpace
   :members:
   :special-members: __call__

.. autoclass:: AdditiveSequence
   :members:
