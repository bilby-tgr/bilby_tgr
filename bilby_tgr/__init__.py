from . import tiger, siqm, mdr, pseob
try:
    from ._version import version as __version__
except ModuleNotFoundError:  # development mode
    __version__ = 'unknown'
