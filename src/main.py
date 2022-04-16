from __future__ import annotations

import os

from underworld import function as fn
from underworld import mpi
from underworld import underworld as uw

from PlatePolygons import SubductionZonePolygons
from modelParameters import (
    ModelParameterBuilder,
    ModelParameterMap,
    ModelParameterDaoBuilder,
    ScalingCoefficient,
    ScalingCoefficientType,
)
from SubductionModel import SubductionModel

u = uw.scaling.units


parameterSet = strak2021ModelParameterSet
# nonDimensionalizeParameters
parameterSet.nonDimensionalizeParameters()

pp = SubductionZonePolygons(
    parameterSet,
    29,
    200e3 * u.meter,
    6000e3 * u.meter,
    30e3 * u.meter,
    20e3 * u.meter,
    30e3 * u.meter,
    100e3 * u.meter,
    500e3 * u.meter,
    200e3 * u.meter,
    3000e3 * u.meter,
    30e3 * u.meter,
    30e3 * u.meter,
    120e3 * u.meter,
)
model = SubductionModel("test", parameterSet, pp)

if __name__ == "__main__":
    outputPath = os.path.join(os.path.abspath("."), "output/")

    if uw.mpi.rank == 0:
        if not os.path.exists(outputPath):
            os.makedirs(outputPath)
    mpi.barrier()
    model.getParticlePlot()
