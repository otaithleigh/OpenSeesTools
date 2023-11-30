# Prefer locally-built openseespy to pip-installed openseespy
try:
    import opensees
except ImportError:
    try:
        from openseespy import opensees
    except ImportError:
        import warnings

        warnings.warn('OpenSeesPy not found on this system.')
        opensees = None

from .basic import *
from .sectionanalysis import *
from .variableanalysis import *
from .wsection import *
