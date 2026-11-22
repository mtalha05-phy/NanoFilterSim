"""NanoFilterSim -- independent Python reference implementation (Transfer Matrix Method)."""
__version__ = "0.1.0"
from .core.tmm import solve, TMMResult
from .materials.database import get_material
