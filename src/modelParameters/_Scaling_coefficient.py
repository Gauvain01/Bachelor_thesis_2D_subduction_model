from abc import ABC, abstractmethod, abstractproperty
from typing import List

from pint.unit import _Unit
from underworld import scaling
from UWGeodynamics import UnitRegistry

from modelParameters._Model_parameter import ModelParameter
from modelParameters._scaling_coefficient_type import ScalingCoefficientType


class ScalingCoefficient(ABC):
    """
    base class for setting and using the Scaling Coefficient algorithm
    """

    def __init__(self) -> None:
        self._setUnderWorldScalingFactors()

    def _setUnderWorldScalingFactors(self):
        sco = scaling.get_coefficients()
        sco["[length]"] = self.lengthCoefficient
        sco["[temperature]"] = self.temperatureCoefficient
        sco["[mass]"] = self.massCoefficient
        sco["[time]"] = self.timeCoefficient

    def nonDimensionalizeUnderworld(self, value):
        return scaling.non_dimensionalise(value)

    @property
    def _functionToCoefficientEnumMapper(self) -> dict:
        mapperDict = {
            ScalingCoefficientType.VELOCITYVECTOR: self.scalingForVelocityVector,
            ScalingCoefficientType.VISCOSITY: self.scalingForViscosity,
            ScalingCoefficientType.STRESS: self.scalingForStress,
            ScalingCoefficientType.TEMPERATURE: self.scalingForTemperature,
            ScalingCoefficientType.LENGTH: self.scalingForLength,
            ScalingCoefficientType.GRADIENT: self.scalingForGradient,
            ScalingCoefficientType.MASS: self.scalingForMass,
            ScalingCoefficientType.TIME: self.scalingForTime,
        }
        return mapperDict

    def nonDimensionalizeUnit(
        self, value: _Unit, scalingTypeEnum: ScalingCoefficientType
    ):
        func = self._functionToCoefficientEnumMapper[scalingTypeEnum]
        return func(value)

    @property
    def viscosityCoefficient(self):
        pass

    @property
    def stressCoefficient(self):
        pass

    @property
    def temperatureCoefficient(self):
        pass

    @property
    def lengthCoefficient(self):
        pass

    @property
    def gradientCoefficient(self):
        pass

    @property
    def velocityVectorCoefficient(self):
        pass

    @property
    def massCoefficient(self):
        pass

    @property
    def timeCoefficient(self):
        pass

    def scalingForViscosity(self, dimensionalValue: UnitRegistry) -> float:
        return (dimensionalValue) / (self.viscosityCoefficient)

    def scalingForStress(self, dimensionalValue: UnitRegistry) -> float:
        return (dimensionalValue) / (self.stressCoefficient)

    def scalingForTemperature(self, dimensionalValue: UnitRegistry) -> float:
        return (dimensionalValue) / (self.temperatureCoefficient)

    def scalingForLength(self, dimensionalValue: UnitRegistry) -> float:
        return (dimensionalValue) / (self.lengthCoefficient)

    def scalingForGradient(self, dimensionalValue: UnitRegistry) -> float:
        return (dimensionalValue) / (self.gradientCoefficient)

    def scalingForVelocityVector(self, dimensionalValue: UnitRegistry) -> float:
        return (dimensionalValue) / (self.velocityVectorCoefficient)

    def scalingForMass(self, dimensionalValue: UnitRegistry) -> float:
        return (dimensionalValue) / (self.massCoefficient)

    def scalingForTime(self, dimensionalValue: UnitRegistry) -> float:
        return (dimensionalValue) / (self.timeCoefficient)
