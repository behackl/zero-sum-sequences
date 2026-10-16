"""Private capability signals, distinct from invalid data or provider failures."""


class _ArithmeticUnavailable(ValueError):
    """A required finite-group capability is absent or explicitly unsupported."""


class _ParentIterationUnavailable(TypeError):
    """The parent explicitly lacks Python iteration protocols."""
