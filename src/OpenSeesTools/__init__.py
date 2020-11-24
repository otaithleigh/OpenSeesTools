# Prefer locally-built openseespy to pip-installed openseespy
try:
    import opensees
except ImportError:
    try:
        from openseespy import opensees
    except ImportError:
        raise RuntimeError('OpenSeesPy not found on this system.')

from .basic import *
from .SectionAnalysis import *
from .VariableAnalysis import *
from .wSection import *
