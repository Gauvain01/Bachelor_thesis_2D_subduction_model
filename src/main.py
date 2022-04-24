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
        name="test2",
        modelParameterMap=get_Strak_2021_model_parameter_map(),
        resolution=(44, 44),
        stepAmountCheckpoint=2,
        subductionZonePolygons=polygons,
        totalSteps=10,
    )
    # model.run()

lv = FigureViewer("./output/test2")
lv.fromCheckPoint(9, model)
