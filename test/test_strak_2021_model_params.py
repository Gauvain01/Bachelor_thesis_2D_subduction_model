from test.test_resources.Strak_2021_model_params_resources import (
    get_Strak_2021_model_parameter_map,
)

import numpy as np
from modelParameters import ModelParameterMapBuilder
from RheologyFunctions import RheologyFunctions
from underworld.scaling import units as u


def test_model_parameter_map():
    parameterMap = get_Strak_2021_model_parameter_map()
    assert parameterMap.gasConstant.nonDimensionalValue.magnitude == 8.3145
    assert parameterMap.modelHeight.nonDimensionalValue.magnitude == 1
    assert parameterMap.referenceDensity.nonDimensionalValue.magnitude == 0
    assert parameterMap.modelLength.nonDimensionalValue.magnitude == 4
    assert (
        format(
            parameterMap.yieldStressOfSpTopLayer.nonDimensionalValue.magnitude, ".2E"
        )
        == f"{5.05e5:.2E}"
    )

    assert parameterMap.thermalDiffusivity.nonDimensionalValue.magnitude == 1e-6
    assert parameterMap.thermalExpansivity.nonDimensionalValue.magnitude == 1e-5
    assert parameterMap.referenceViscosity.nonDimensionalValue.magnitude == 1
    assert parameterMap.referenceTemperature.nonDimensionalValue.magnitude == 1
    assert parameterMap.upperMantleViscosity.nonDimensionalValue.magnitude == 1
    assert parameterMap.lowerMantleViscosity.nonDimensionalValue.magnitude == 100
    assert parameterMap.spBottomLayerViscosity.nonDimensionalValue.magnitude == 50
    assert (
        round(parameterMap.spCoreLayerViscosity.nonDimensionalValue.magnitude) == 1000
    )
    assert round(parameterMap.spTopLayerViscosity.nonDimensionalValue.magnitude) == 1000


def test_blueprint_map():
    bluePrint = get_Strak_2021_model_parameter_map(True)

    builder = ModelParameterMapBuilder.fromBluePrint(bluePrint)
    paramMap = builder.build()
    assert repr(paramMap) == repr(get_Strak_2021_model_parameter_map())


def test_blueprint_map_with_change():
    bluePrint = get_Strak_2021_model_parameter_map(True)

    builder = ModelParameterMapBuilder.fromBluePrint(bluePrint)
    builder.setGasConstant(u.Quantity(100.0))
    paramMap = builder.build()
    assert repr(paramMap) != repr(get_Strak_2021_model_parameter_map())
    assert paramMap.gasConstant.nonDimensionalValue.magnitude == 100.0


def test_rayleigh_number():
    rheologyFn = RheologyFunctions(get_Strak_2021_model_parameter_map(), None)
    rayleighNumber = rheologyFn.getRayleighNumber()
    print(rayleighNumber)
    assert format(rayleighNumber, ".1E") == format(3.5e7, ".1E")
