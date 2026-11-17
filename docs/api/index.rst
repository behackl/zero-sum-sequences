API reference
=============

.. py:module:: zero_sum_sequences

Everything is imported from :mod:`zero_sum_sequences`. A typical session
builds a :class:`FiniteAdditiveGroup` and an :class:`AdditiveSequenceSpace`
over it, constructs sequences, and calls their methods, such as
:meth:`~AdditiveSequence.length_set` or :meth:`~AdditiveSequence.orbit`.
The remaining classes are for reusing work across many queries and for
finer control.

.. toctree::
   :maxdepth: 1

   sequences
   factorizations
   symmetry
   constants
