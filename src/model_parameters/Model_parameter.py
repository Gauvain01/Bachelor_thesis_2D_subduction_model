from __future__ import annotations

from dataclasses import dataclass
from typing import Union

from pint.quantity import _Quantity
from pint.unit import _Unit

from model_parameters.scaling_coefficient_type_enum import ScalingCoefficientTypeEnum


@dataclass(frozen=True, repr=True)
class ModelParameter:
    dimensionalValue: Union[_Unit, _Quantity]
    nonDimensionalValue: float
