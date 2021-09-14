#!/usr/bin/env/python
"""Collection of data structure utilities."""

import collections
import sys

# https://stackoverflow.com/questions/4126348/how-can-this-function-be-rewritten-to-implement-ordereddict/4127426#4127426
# Thanks to martineau

class OrderedDefaultdict(collections.OrderedDict):
    """A defaultdict with OrderedDict as its base class."""

    def __init__(self, default_factory=None, *args, **kwargs):
        if not (default_factory is None
                or isinstance(default_factory, collections.Callable)):
            raise TypeError('first argument must be callable or None')
        super(OrderedDefaultdict, self).__init__(*args, **kwargs)
        self.default_factory = default_factory  # called by __missing__()

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key, )
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):  # optional, for pickle support
        args = (self.default_factory,) if self.default_factory else tuple()
        return self.__class__, args, None, None, self.items()

    def __repr__(self):  # optional
        return '%s(%r, %r)' % (self.__class__.__name__, self.default_factory,
                               list(self.items()))

