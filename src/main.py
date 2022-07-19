from __future__ import annotations

from underworld import visualisation
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
        3000e3 * u.meter,
        30e3 * u.meter,
        20e3 * u.meter,
        30e3 * u.meter,
        100e3 * u.meter,
    )
    print(
        get_Strak_2021_model_parameter_map().yieldStressOfSpTopLayer.nonDimensionalValue.magnitude
    )
    model = SubductionModel(
        modelParameters=get_Strak_2021_model_parameter_map(),
        totalSteps=10,
        checkPointSteps=2,
        resolution=(300, 100),
        subductionZonePolygons=polygons,
        name="test_34",
    )
    print("banaan")
    model.run()

    # # fig = visualisation.Figure(figsize=(3000, 1000))
    # # fig.append(visualisation.objects.Mesh(model.mesh))
    # # fig.show()
    # lv = FigureViewer("src/output/test_16")
    # lv.fromDb(100)
