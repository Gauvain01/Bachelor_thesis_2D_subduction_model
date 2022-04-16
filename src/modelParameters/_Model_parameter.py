from __future__ import annotations

import attr
from pint.quantity import _Quantity


@attr.s(frozen=True, repr=True, slots=True, kw_only=True)
class ModelParameter:
    dimensionalValue: _Quantity = attr.ib(
        validator=attr.validators.instance_of(_Quantity)
    )
    nonDimensionalValue: _Quantity = attr.ib(
        validator=attr.validators.instance_of(_Quantity)
    )
