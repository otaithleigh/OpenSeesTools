# Prefer locally-built openseespy to pip-installed openseespy
try:
    import opensees
except ImportError:
    from openseespy import opensees

from .basic import *
from .SectionAnalysis import *
from .VariableAnalysis import *
from .wSection import *
