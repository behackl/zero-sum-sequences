"""Tools for finite additive sequences and zero-sum computations."""

from importlib.metadata import version as _distribution_version

from .additive_sequence import AdditiveSequence, AdditiveSequenceSpace
from .atom_catalogue import AtomCatalogue, AtomOrbit
from .automorphisms import (
    Automorphism,
    AutomorphismGroup,
    AutomorphismGroupUnavailable,
    CanonicalForm,
)
from .factorization import Factorization, FactorizationSolver, FactorizationStatistics
from .factorization_cache import (
    AnchoredCanonicalizer,
    FactorizationCache,
    FactorizationCacheStatistics,
)
from .invariants import (
    ZeroSumConstant,
    davenport_constant,
    erdos_ginzburg_ziv_constant,
    eta_constant,
    exponent,
    zero_sum_constant,
)
from .orbits import (
    AutomorphismAction,
    AutomorphismActionUnavailable,
    OrbitWitness,
)
from .parents import FiniteAdditiveGroup
from .relations import FactorizationRelation

__all__ = [
    "AdditiveSequence",
    "AdditiveSequenceSpace",
    "AtomCatalogue",
    "AtomOrbit",
    "Factorization",
    "FactorizationRelation",
    "AnchoredCanonicalizer",
    "FactorizationCache",
    "FactorizationCacheStatistics",
    "FactorizationSolver",
    "FactorizationStatistics",
    "FiniteAdditiveGroup",
    "Automorphism",
    "AutomorphismAction",
    "AutomorphismActionUnavailable",
    "AutomorphismGroup",
    "AutomorphismGroupUnavailable",
    "CanonicalForm",
    "OrbitWitness",
    "ZeroSumConstant",
    "davenport_constant",
    "erdos_ginzburg_ziv_constant",
    "eta_constant",
    "exponent",
    "zero_sum_constant",
]
__version__ = _distribution_version("zero-sum-sequences")
