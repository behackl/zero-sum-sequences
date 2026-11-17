Automorphisms and orbits
========================

.. currentmodule:: zero_sum_sequences

Automorphisms of the group act on sequences term by term. For groups from
:meth:`FiniteAdditiveGroup.cyclic_product`, the full automorphism group is
available without setup, and :meth:`AdditiveSequence.orbit` and
:meth:`AdditiveSequenceSpace.automorphism_group` are the usual entry
points. The classes below are needed for custom groups of automorphisms
and for canonical forms.

.. autosummary::
   :nosignatures:

   AutomorphismGroup
   Automorphism
   CanonicalForm
   OrbitWitness
   AutomorphismAction

.. autoclass:: AutomorphismGroup
   :members:

.. autoclass:: Automorphism
   :members:

.. autoclass:: CanonicalForm
   :members:

.. autoclass:: OrbitWitness
   :members:

.. autoclass:: AutomorphismAction
   :members:

.. autoexception:: AutomorphismGroupUnavailable

.. autoexception:: AutomorphismActionUnavailable
