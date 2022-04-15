from __future__ import annotations

from typing import Union

from pint.quantity import _Quantity
from pint.unit import _Unit

from modelParameters._Model_parameter import ModelParameter
from modelParameters._Model_parameter_builder import ModelParameterBuilder
from modelParameters._Model_parameter_dao import ModelParameterDao
from modelParameters._scaling_coefficient_type import ScalingCoefficientType


class ModelParameterDaoBuilder:
    def __init__(self, modelParameterBuilder: ModelParameterBuilder) -> None:
        self._modelParameterbuilder = modelParameterBuilder
        self._modelHeight: ModelParameter
        self._modelLength: ModelParameter
        self._referenceDensity: ModelParameter
        self._gravitationalAcceleration: ModelParameter
        self._referenceTemperature: ModelParameter
        self._thermalExpansivity: ModelParameter
        self._thermalDiffusivity: ModelParameter
        self._referenceViscosity: ModelParameter
        self._upperMantleViscosity: ModelParameter
        self._lowerMantleViscosity: ModelParameter
        self._spTopLayerViscosity: ModelParameter
        self._spCoreLayerViscosity: ModelParameter
        self._spBottomLayerViscosity: ModelParameter
        self._yieldStressOfSpTopLayer: ModelParameter
        self._gasConstant: ModelParameter
        self._lowerMantleHeigth: ModelParameter
        self._coreShearModulus: ModelParameter
        self._timeScaleStress: ModelParameter

    def _reset(self):
        self._modelHeight: ModelParameter = None
        self._modelLength: ModelParameter = None
        self._referenceDensity: ModelParameter = None
        self._gravitationalAcceleration: ModelParameter = None
        self._referenceTemperature: ModelParameter = None
        self._thermalExpansivity: ModelParameter = None
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
        self._coreShearModulus: ModelParameter = None
        self._timeScaleStress: ModelParameter = None

    def checkAttribute(self, attribute):
        if attribute is not None:
            raise AttributeError(f"{attribute = } already assigned")

    def setModelHeight(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
    ) -> ModelParameterDaoBuilder:
        self._checkAttribute(self._modelHeight)
        self._modelHeight = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum
        )
        return self

    def setModelLength(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
    ) -> ModelParameterDaoBuilder:
        self._checkAttribute(self._modelLength)
        self._modelLength = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum
        )
        return self

    def setReferenceDensity(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
    ) -> ModelParameterDaoBuilder:
        self._checkAttribute(self._referenceDensity)
        self._referenceDensity = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum
        )
        return self

    def setGravitationalAcceleration(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
    ) -> ModelParameterDaoBuilder:
        self._checkAttribute(self._gravitationalAcceleration)
        self._gravitationalAcceleration = (
            self._modelParameterbuilder.buildModelParameter(
                value, scalingCoefficientTypeEnum
            )
        )
        return self

    def setReferenceTemperature(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
    ) -> ModelParameterDaoBuilder:
        self._checkAttribute(self._referenceTemperature)
        self._referenceTemperature = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum
        )
        return self

    def setThermalExpansivity(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
    ) -> ModelParameterDaoBuilder:
        self._checkAttribute(self._thermalExpansivity)
        self._thermalExpansivity = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum
        )
        return self

    def setThermalDiffusivity(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
    ) -> ModelParameterDaoBuilder:
        self._checkAttribute(self._thermalDiffusivity)
        self._thermalDiffusivity = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum
        )
        return self

    def setReferenceViscosity(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
    ) -> ModelParameterDaoBuilder:
        self._checkAttribute(self._referenceViscosity)
        self._referenceViscosity = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum
        )
        return self

    def setUpperMantleViscosity(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
    ) -> ModelParameterDaoBuilder:
        self._checkAttribute(self._upperMantleViscosity)
        self._upperMantleViscosity = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum
        )
        return self

    def setLowerMantleViscosity(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
    ) -> ModelParameterDaoBuilder:
        self._checkAttribute(self._lowerMantleViscosity)
        self._lowerMantleViscosity = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum
        )
        return self

    def setSpTopLayerViscosity(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
    ) -> ModelParameterDaoBuilder:
        self._checkAttribute(self._spTopLayerViscosity)
        self._spTopLayerViscosity = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum
        )
        return self

    def setSpCoreLayerViscosity(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
    ) -> ModelParameterDaoBuilder:
        self._checkAttribute(self._spCoreLayerViscosity)
        self._spCoreLayerViscosity = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum
        )
        return self

    def setSpBottomLayerViscosity(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
    ) -> ModelParameterDaoBuilder:
        self._checkAttribute(self._spBottomLayerViscosity)
        self._spBottomLayerViscosity = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum
        )
        return self

    def setYieldStressOfSpTopLayer(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
    ) -> ModelParameterDaoBuilder:
        self._checkAttribute(self._yieldStressOfSpTopLayer)
        self._yieldStressOfSpTopLayer = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum
        )
        return self

    def setGasConstant(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
    ) -> ModelParameterDaoBuilder:
        self._checkAttribute(self._gasConstant)
        self._gasConstant = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum
        )
        return self

    def setLowerMantleHeight(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
    ) -> ModelParameterDaoBuilder:
        self._checkAttribute(self._lowerMantleHeigth)
        self._lowerMantleHeigth = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum
        )
        return self

    def setCoreShearModulus(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
    ) -> ModelParameterDaoBuilder:
        self._checkAttribute(self._coreShearModulus)
        self._coreShearModulus = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum
        )
        return self

    def setTimeScaleStress(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientType] = None,
    ) -> ModelParameterDaoBuilder:
        self._checkAttribute(self._timeScaleStress)
        self._timeScaleStress = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum
        )
        return self

    def build(self) -> ModelParameterDao:
        modelParameterDao = ModelParameterDao(
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
            timeScaleStress=self._timeScaleStress,
        )
        self._reset()
        return modelParameterDao
