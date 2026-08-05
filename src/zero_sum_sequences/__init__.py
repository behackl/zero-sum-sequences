"""Tools for finite additive sequences and zero-sum computations."""

from .additive_sequence import AdditiveSequence, AdditiveSequenceSpace
from .atom_catalogue import AtomCatalogue
from .factorization import FactorizationSolver

__all__ = [
    "AdditiveSequence",
    "AdditiveSequenceSpace",
    "AtomCatalogue",
    "FactorizationSolver",
]
__version__ = "0.1.0"
