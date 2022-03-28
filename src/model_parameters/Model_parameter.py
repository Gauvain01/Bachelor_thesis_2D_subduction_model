from __future__ import annotations

import UWGeodynamics as geo
from UWGeodynamics import UnitRegistry

from model_parameters.scaling_coefficient_type_enum import ScalingCoefficientTypeEnum


class ModelParameter:
    def __init__(
        self,
        name: str,
        dimensionalValue: UnitRegistry,
        scalingCoefficientTypeEnum: ScalingCoefficientTypeEnum,
    ) -> None:
        self.name = name
        self.dimensionalValue = dimensionalValue
        self.scalingCoefficientTypeEnum = scalingCoefficientTypeEnum
        self.nonDimensionalValue: float = None

    def setNondimensionalValue(self, value: float) -> ModelParameter:
        self.nonDimensionalValue = value
        return self
