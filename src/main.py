from __future__ import annotations

from underworld.scaling import units as u

from FigureViewer import FigureViewer
from PlatePolygons import SubductionZonePolygons
from strakParam import get_Strak_2021_model_parameter_map
from OldSubductionModel import OldSubductionModel

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
    model = OldSubductionModel(
        name="test2",
        modelParameterMap=get_Strak_2021_model_parameter_map(),
        resolution=(200, 100),
        stepAmountCheckpoint=50,
        subductionZonePolygons=polygons,
        totalSteps=100,
    )
    # model.run()

lv = FigureViewer("./output/test2")
lv.fromCheckPoint(50, model)
