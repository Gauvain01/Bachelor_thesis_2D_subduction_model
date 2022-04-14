from typing import Union

from pint.quantity import _Quantity
from pint.unit import _Unit

from modelParameters._Model_parameter import ModelParameter
from modelParameters._Scaling_coefficient import ScalingCoefficient
from modelParameters._scaling_coefficient_type import ScalingCoefficientType


class ModelParameterBuilder:
    def __init__(self, scalingCoefficient: ScalingCoefficient) -> None:
        self.scalingCoefficient = scalingCoefficient

    def buildModelParameter(
        self,
        value: Union[_Unit, _Quantity],
        scalingCoefficientEnum: Union[None, ScalingCoefficientType] = None,
    ) -> ModelParameter:
        if isinstance(value, _Unit):
            if scalingCoefficientEnum is None:
                raise ValueError(
                    "value type is pint.unit._Unit but no scalingCoefficientType Was Given"
                )
            elif scalingCoefficientEnum == ScalingCoefficientType.NONE:
                nonDimensionalValue = value.magnitude
            else:
                nonDimensionalValue = self.scalingCoefficient.nonDimensionalizeUnit(
                    value, scalingCoefficientEnum
                )
        else:
            nonDimensionalValue = value.magnitude
        return ModelParameter(value, nonDimensionalValue)
