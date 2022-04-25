from __future__ import annotations

import inspect
from copy import deepcopy
from os import name
from typing import Any, Dict, Union

from attr import attributes
from pint.quantity import _Quantity
from pint.unit import _Unit
from typing_extensions import Self

from modelParameters._Model_parameter import ModelParameter
from modelParameters._Model_parameter_builder import ModelParameterBuilder
from modelParameters._Model_parameter_map import ModelParameterMap
from modelParameters._scaling_coefficient_type import ScalingCoefficientType


class PassableException(Exception):
    pass


class AlreadyAssignedError(Exception):
    pass


class ModelParameterMapBuilder:
    def __init__(self, modelParameterBuilder: ModelParameterBuilder) -> None:
        self._modelParameterbuilder = modelParameterBuilder
        self._modelHeight: ModelParameter = None
        self._modelLength: ModelParameter = None
        self._referenceDensity: ModelParameter = None
        self._gravitationalAcceleration: ModelParameter = None
        self._referenceTemperature: ModelParameter = None
        self._defaultStrainRate: ModelParameter = None
        self._thermalExpansivity: ModelParameter = None
        self._minimalStrainRate: ModelParameter = None
        self._thermalDiffusivity: ModelParameter = None
        self._referenceViscosity: ModelParameter = None
        self._upperMantleViscosity: ModelParameter = None
        self._lowerMantleViscosity: ModelParameter = None
        self._spTopLayerViscosity: ModelParameter = None
        self._spCoreLayerViscosity: ModelParameter = None
        self._spBottomLayerViscosity: ModelParameter = None
        self._yieldStressOfSpTopLayer: ModelParameter = None
        self._gasConstant: ModelParameter = None
        self._lowerMantleHeigth: ModelParameter = None
        self._coreShearModulus = None
        self._timeScaleStress: ModelParameter = None
        self._fromBluePrint = False
        self._temperatureContrast = None
        self._resetFlag = False

    def _reset(self):
        self._resetFlag = True
        try:
            for name, value in inspect.getmembers(self):
                if not name.startswith("__") and not inspect.ismethod(value):
                    if isinstance(value, ModelParameter) or name == "_bluePrint":
                        self.__setattr__(name, None)
        finally:
            self._resetFlag = False

    def __setattr__(self, __name: str, __value: Any) -> None:
        try:
            if not self._resetFlag:
                attribute = self.__getattribute__(__name)
                if isinstance(attribute, ModelParameter):
                    raise AlreadyAssignedError(
                        f"{__name}{attribute = } already assigned"
                    )
                else:
                    super().__setattr__(__name, __value)
            else:
                super().__setattr__(__name, __value)
        except AttributeError:
            super().__setattr__(__name, __value)

    def setTemperatureContrast(
        self,
        value: _Quantity,
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType],
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._temperatureContrast = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
        )
        return self

    def setModelHeight(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._modelHeight = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
        )
        return self

    def setMinimalStrainRate(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._minimalStrainRate = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
        )
        return self

    def setModelLength(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._modelLength = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
        )
        return self

    def setReferenceDensity(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._referenceDensity = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
        )
        return self

    def setGravitationalAcceleration(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._gravitationalAcceleration = (
            self._modelParameterbuilder.buildModelParameter(
                value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
            )
        )
        return self

    def setReferenceTemperature(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._referenceTemperature = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
        )
        return self

    def setThermalExpansivity(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._thermalExpansivity = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
        )
        return self

    def setThermalDiffusivity(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._thermalDiffusivity = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
        )
        return self

    def setReferenceViscosity(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._referenceViscosity = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
        )
        return self

    def setUpperMantleViscosity(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._upperMantleViscosity = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
        )
        return self

    def setLowerMantleViscosity(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._lowerMantleViscosity = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
        )
        return self

    def setSpTopLayerViscosity(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._spTopLayerViscosity = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
        )
        return self

    def setSpCoreLayerViscosity(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._spCoreLayerViscosity = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
        )
        return self

    def setSpBottomLayerViscosity(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._spBottomLayerViscosity = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
        )
        return self

    def setYieldStressOfSpTopLayer(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._yieldStressOfSpTopLayer = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
        )
        return self

    def setGasConstant(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._gasConstant = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
        )
        return self

    def setLowerMantleHeight(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._lowerMantleHeigth = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
        )
        return self

    def setDefaultStrainRate(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._defaultStrainRate = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
        )
        return self

    def setCoreShearModulus(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._coreShearModulus = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
        )
        return self

    def setDeltaTime(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
        nonDimensionalOverrideValue=None,
    ) -> ModelParameterMapBuilder:
        self._deltaTime = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum, nonDimensionalOverrideValue
        )
        return self

    def createBluePrint(self):
        output = {}
        for name, value in inspect.getmembers(self):
            if not name.startswith("__") and not inspect.ismethod(value):
                if isinstance(value, ModelParameter):
                    output[name] = value
        output["modelParameterBuilder"] = self._modelParameterbuilder
        return output

    @classmethod
    def fromBluePrint(cls, bluePrint: Dict) -> ModelParameterMapBuilder:
        parameterBuilder = bluePrint.pop("modelParameterBuilder")
        obj = cls(parameterBuilder)
        obj._fromBluePrint = True
        obj._bluePrint = bluePrint
        return obj

    def build(self) -> ModelParameterMap:
        if self._fromBluePrint is False:
            modelParameterDao = ModelParameterMap(
                scalingCoefficient=self._modelParameterbuilder.scalingCoefficient,
                modelHeight=self._modelHeight,
                modelLength=self._modelLength,
                referenceDensity=self._referenceDensity,
                gravitationalAcceleration=self._gravitationalAcceleration,
                referenceTemperature=self._referenceTemperature,
                thermalExpansivity=self._thermalExpansivity,
                thermalDiffusivity=self._thermalDiffusivity,
                referenceViscosity=self._referenceViscosity,
                upperMantleViscosity=self._upperMantleViscosity,
                lowerMantleViscosity=self._lowerMantleViscosity,
                spTopLayerViscosity=self._spTopLayerViscosity,
                spCoreLayerViscosity=self._spCoreLayerViscosity,
                spBottomLayerViscosity=self._spBottomLayerViscosity,
                yieldStressOfSpTopLayer=self._yieldStressOfSpTopLayer,
                gasConstant=self._gasConstant,
                lowerMantleHeigth=self._lowerMantleHeigth,
                deltaTime=self._deltaTime,
                coreShearModulus=self._coreShearModulus,
                temperatureContrast=self._temperatureContrast,
                defaultStrainRate=self._defaultStrainRate,
                minimalStrainRate=self._minimalStrainRate,
            )
            self._reset()
            return modelParameterDao
        else:
            for name, value in self._bluePrint.items():
                obj = self.__getattribute__(name)

                if obj is None:
                    self.__setattr__(name, value)

            self._fromBluePrint = False
            return self.build()
