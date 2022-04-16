from test.test_resources.Strak_2021_model_params_resources import (
    get_Strak_2021_model_parameter_map,
)

import numpy as np
from modelParameters import ModelParameterMapBuilder
from RheologyFunctions import RheologyFunctions


def test_model_parameter_map():
    strakDao = get_Strak_2021_model_parameter_map()
    assert strakDao.gasConstant.nonDimensionalValue == 8.3145
    assert strakDao.modelHeight.nonDimensionalValue == 1
    assert strakDao.referenceDensity.nonDimensionalValue == 0
    assert strakDao.modelLength.nonDimensionalValue == 4
    assert (
        format(strakDao.yieldStressOfSpTopLayer.nonDimensionalValue, ".2E")
        == f"{5.05e5:.2E}"
    )

    assert strakDao.thermalDiffusivity.nonDimensionalValue == 1e-6
    assert strakDao.thermalExpansivity.nonDimensionalValue == 1e-5
    assert strakDao.referenceViscosity.nonDimensionalValue == 1
    assert strakDao.referenceTemperature.nonDimensionalValue == 1
    assert strakDao.upperMantleViscosity.nonDimensionalValue == 1
    assert strakDao.lowerMantleViscosity.nonDimensionalValue == 100
    assert strakDao.spBottomLayerViscosity.nonDimensionalValue == 50
    assert round(strakDao.spCoreLayerViscosity.nonDimensionalValue) == 1000
    assert round(strakDao.spTopLayerViscosity.nonDimensionalValue) == 1000


def test_blueprint_dao():
    bluePrint = get_Strak_2021_model_parameter_map(True)

    builder = ModelParameterMapBuilder.fromBluePrint(bluePrint)
    dao = builder.build()
    assert repr(dao) == repr(get_Strak_2021_model_parameter_map())


def test_rayleigh_number():
    rheologyFn = RheologyFunctions(get_Strak_2021_model_parameter_map(), None)
    rayleighNumber = rheologyFn.getRayleighNumber()
    print(rayleighNumber)
    assert format(rayleighNumber, ".1E") == format(3.5e7, ".1E")
