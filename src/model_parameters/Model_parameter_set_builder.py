from __future__ import annotations
from typing import Union
from model_parameters.Model_parameter_builder import ModelParameterBuilder
from model_parameters.Model_parameter import ModelParameter
from model_parameters.scaling_coefficient_type_enum import ScalingCoefficientTypeEnum
from pint.unit import _Unit
from pint.quantity import _Quantity


class ModelParameterSetBuilder:
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

    def setModelHeight(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientTypeEnum] = None,
    ) -> ModelParameterSetBuilder:
        self._modelHeight = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum
        )
        return self

    def setModelLength(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientTypeEnum: Union[None, ScalingCoefficientTypeEnum] = None,
    ) -> ModelParameterSetBuilder:
        self._modelLength = self._modelParameterbuilder.buildModelParameter(
            value, scalingCoefficientTypeEnum
        )
        return self
