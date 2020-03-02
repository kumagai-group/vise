# -*- coding: utf-8 -*-
from custodian.vasp import handlers as orig_handlers

from vise.custodian_extension.error_handlers import (
    ViseVaspErrorHandler, ViseUnconvergedErrorHandler, MemoryOverflowHandler,
    DivergingEnergyErrorHandler, TooLongTimeCalcErrorHandler)
from vise.util.logger import get_logger

# Note1: Custodian handler groups are modified to avoid IBRION=1, which does not
#        show efermi (& eigenvalues?) in vasprun.xml.
#        Therefore, IBRION=1 could be switched on only for rough calculations.
# Note2: Don't forget "()" to generate the ErrorHandler class object.

logger = get_logger(__name__)


def handler_group(name: str = "default",
                  nionic_steps: int = 1,
                  timeout: int = 518400) -> list:

    handlers = {
        "minimum":    [orig_handlers.MeshSymmetryErrorHandler(),
                       orig_handlers.PotimErrorHandler(),
                       orig_handlers.PositiveEnergyErrorHandler(),
                       orig_handlers.FrozenJobErrorHandler(),
                       orig_handlers.StdErrHandler(),
                       ViseVaspErrorHandler(),
                       ViseUnconvergedErrorHandler(),
                       MemoryOverflowHandler()],
        "default":    [orig_handlers.MeshSymmetryErrorHandler(),
                       orig_handlers.NonConvergingErrorHandler(),
                       orig_handlers.PotimErrorHandler(),
                       orig_handlers.PositiveEnergyErrorHandler(),
                       orig_handlers.FrozenJobErrorHandler(),
                       orig_handlers.StdErrHandler(),
                       ViseVaspErrorHandler(),
                       ViseUnconvergedErrorHandler(),
                       MemoryOverflowHandler(),
                       DivergingEnergyErrorHandler(),
                       TooLongTimeCalcErrorHandler(timeout=timeout),
                       ],
        "dielectric": [orig_handlers.MeshSymmetryErrorHandler(),
                       orig_handlers.NonConvergingErrorHandler(
                           nionic_steps=nionic_steps),
                       orig_handlers.PositiveEnergyErrorHandler(),
                       orig_handlers.FrozenJobErrorHandler(),
                       orig_handlers.StdErrHandler(),
                       orig_handlers.LrfCommutatorHandler(),
                       ViseUnconvergedErrorHandler(),
                       MemoryOverflowHandler(),
                       DivergingEnergyErrorHandler(),
                       TooLongTimeCalcErrorHandler(timeout=timeout),
                       ],
        "no_handler": []
    }

    try:
        return handlers[name]
    except KeyError:
        logger.critical(f"Accepted handler group names are  {handlers.keys()}")
        raise

