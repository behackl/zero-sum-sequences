"""Pickle support that preserves the identity of configured parents and spaces.

Sequences compare only within one space object. A pickled group or space
therefore carries a token: unpickling in a process that still holds the
original (or an earlier copy) returns that object, so results sent back from
worker processes belong to the caller's space.
"""

from __future__ import annotations

import uuid
import weakref

_live: weakref.WeakValueDictionary[str, object] = weakref.WeakValueDictionary()


def _token(obj) -> str:
    """Return the object's pickle token, registering it on first use."""

    token = obj._pickle_token
    if token is None:
        token = uuid.uuid4().hex
        obj._pickle_token = token
        _live[token] = obj
    return token


def _restore(token: str, build, args: tuple, kwargs: dict, restore_state=None, state=None):
    """Return the live object for ``token``, or build and register a new one."""

    obj = _live.get(token)
    if obj is None:
        obj = build(*args, **kwargs)
        obj._pickle_token = token
        _live[token] = obj
        if restore_state is not None:
            restore_state(obj, state)
    return obj
