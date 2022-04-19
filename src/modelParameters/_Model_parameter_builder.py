from typing import Union

from pint.quantity import _Quantity
from pint.unit import _Unit
from underworld.scaling import units as u

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
        nonDimensionalOverrideValue=None,
    ) -> ModelParameter:

        if isinstance(value, _Quantity):
            if scalingCoefficientEnum is None and nonDimensionalOverrideValue is None:
                nonDimensionalValue = u.Quantity(value.magnitude)
            elif scalingCoefficientEnum == ScalingCoefficientType.NONE:
                nonDimensionalValue = u.Quantity(value.magnitude)
            elif nonDimensionalOverrideValue is not None:
                if not isinstance(
                    nonDimensionalOverrideValue, float
                ) and not isinstance(nonDimensionalOverrideValue, _Quantity):
                    raise ValueError
                elif isinstance(nonDimensionalOverrideValue, _Quantity):
                    nonDimensionalValue = nonDimensionalOverrideValue
                else:
                    nonDimensionalValue = u.Quantity(nonDimensionalOverrideValue)
            else:
                nonDimensionalValue = self.scalingCoefficient.nonDimensionalizeUnit(
                    value, scalingCoefficientEnum
                )
            return ModelParameter(
                dimensionalValue=value,
                nonDimensionalValue=nonDimensionalValue,
            )
        else:
            raise ValueError
