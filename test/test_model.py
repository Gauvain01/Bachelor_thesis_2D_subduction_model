from test.test_resources.Strak_2021_model_params_resources import (
    get_Strak_2021_model_parameter_map,
)

from PlatePolygons import SubductionZonePolygons
from SubductionModel import SubductionModel
from underworld.scaling import units as u


def test_model():
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
        name="test",
        modelParameterMap=get_Strak_2021_model_parameter_map(),
        stepAmountCheckpoint=10,
        subductionZonePolygons=polygons,
        totalSteps=100,
    )

    assert isinstance(polygons, SubductionZonePolygons)
    assert isinstance(model, SubductionModel)
