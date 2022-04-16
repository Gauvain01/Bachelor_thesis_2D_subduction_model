from __future__ import annotations

import attr
from typing import Union
from pint.quantity import _Quantity
from pint.unit import _Unit


@attr.s(frozen=True, repr=True, slots=True, kw_only=True)
class ModelParameter:
    dimensionalValue: Union[_Unit, _Quantity] = attr.ib(
        validator=attr.validators.instance_of((_Unit, _Quantity))
    )
    nonDimensionalValue: float = attr.ib(validator=attr.validators.instance_of(float))
