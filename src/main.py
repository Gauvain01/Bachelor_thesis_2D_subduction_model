from __future__ import annotations

from underworld.scaling import units as u

from FigureViewer import FigureViewer
from PlatePolygons import SubductionZonePolygons
from strakParam import get_Strak_2021_model_parameter_map
from SubductionModel import SubductionModel

if __name__ == "__main__":
    polygons = SubductionZonePolygons(
        get_Strak_2021_model_parameter_map(),
        27,
        200e3 * u.meter,
        6000e3 * u.meter,
        30e3 * u.meter,
        20e3 * u.meter,
        30e3 * u.meter,
        100e3 * u.meter,
    )
    model = SubductionModel(
        modelParameters=get_Strak_2021_model_parameter_map(),
        totalSteps=10,
        checkPointSteps=2,
        resolution=(40, 40),
        subductionZonePolygons=polygons,
        name="test",
    )
    model.run()

    # lv = FigureViewer("./output/test")
    # lv.fromDb(0)
