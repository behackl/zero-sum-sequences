Automorphisms and orbits
========================

.. currentmodule:: zero_sum_sequences

Automorphisms act on sequences term by term. Groups from
:meth:`FiniteAdditiveGroup.cyclic_product` provide their full automorphism
group, available as :meth:`AdditiveSequenceSpace.automorphism_group`;
:meth:`AdditiveSequence.orbit` and related methods use it by default.

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
