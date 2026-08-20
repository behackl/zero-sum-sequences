"""Tools for finite additive sequences and zero-sum computations."""

from importlib.metadata import version as _distribution_version

from .additive_sequence import AdditiveSequence, AdditiveSequenceSpace
from .atom_catalogue import AtomCatalogue
from .factorization import FactorizationSolver
from .orbits import (
    AutomorphismAction,
    AutomorphismActionUnavailable,
    OrbitWitness,
)
from .parents import FiniteAdditiveGroup

__all__ = [
    "AdditiveSequence",
    "AdditiveSequenceSpace",
    "AtomCatalogue",
    "FactorizationSolver",
    "FiniteAdditiveGroup",
    "AutomorphismAction",
    "AutomorphismActionUnavailable",
    "OrbitWitness",
]
__version__ = _distribution_version("zero-sum-sequences")
