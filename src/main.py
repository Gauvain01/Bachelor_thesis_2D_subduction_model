from __future__ import annotations

from underworld import mpi, visualisation
from underworld.scaling import units as u

from DipDataCollector import DipDataCollector
from FigureViewer import FigureViewer
from PlatePolygons import SubductionZonePolygons
from strakParam import get_Strak_2021_model_parameter_map
from SubductionModel import SubductionModel
from TracerParticle import TracerParticle

if __name__ == "__main__":

    polygons = SubductionZonePolygons(
        get_Strak_2021_model_parameter_map(),
        27,
        200e3 * u.meter,
        3000e3 * u.meter,
        30e3 * u.meter,
        20e3 * u.meter,
        30e3 * u.meter,
        100e3 * u.meter,
    )

    model = SubductionModel(
        modelParameters=get_Strak_2021_model_parameter_map(),
        totalSteps=200,
        checkPointSteps=20,
        resolution=(300, 100),
        subductionZonePolygons=polygons,
        name=f"test_94",
    )
    model.addTracer(TracerParticle(polygons.getHingeCoordinate(), "hinge"))
    model.addTracer(
        TracerParticle(polygons.getLowestPointCoreSlabCoord(), "lowest_point")
    )
    model.addTracer(
        TracerParticle(polygons.getMiddlePartUpperSlabCoord(), "middle_part_upper_slab")
    )
    model.addDipDataCollector(
        DipDataCollector(0.0, 100e3, get_Strak_2021_model_parameter_map(), "shallow")
    )
    model.addDipDataCollector(
        DipDataCollector(150e3, 400e3, get_Strak_2021_model_parameter_map(), "mid")
    )

    model.run()
