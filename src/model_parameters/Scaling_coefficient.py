from abc import ABC, abstractmethod, abstractproperty
from typing import List

from underworld import scaling
from UWGeodynamics import UnitRegistry

from model_parameters.Model_parameter import ModelParameter
from model_parameters.scaling_coefficient_type_enum import ScalingCoefficientTypeEnum


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
            ScalingCoefficientTypeEnum.VELOCITYVECTOR: self.scalingForVelocityVector,
            ScalingCoefficientTypeEnum.VISCOSITY: self.scalingForViscosity,
            ScalingCoefficientTypeEnum.STRESS: self.scalingForStress,
            ScalingCoefficientTypeEnum.TEMPERATURE: self.scalingForTemperature,
            ScalingCoefficientTypeEnum.LENGTH: self.scalingForLength,
            ScalingCoefficientTypeEnum.GRADIENT: self.scalingForGradient,
            ScalingCoefficientTypeEnum.MASS: self.scalingForMass,
            ScalingCoefficientTypeEnum.TIME: self.scalingForTime,
        }
        return mapperDict

    def nonDimensionalizeModelParameter(self, modelParameter: ModelParameter) -> None:

        if (
            modelParameter.scalingCoefficientTypeEnum
            == ScalingCoefficientTypeEnum.MANUAL
        ):
            if modelParameter.nonDimensionalValue == None:
                raise ValueError
        elif (
            modelParameter.scalingCoefficientTypeEnum == ScalingCoefficientTypeEnum.NONE
        ):
            if modelParameter.nonDimensionalValue != None:
                raise ValueError
        else:
            func = self._functionToCoefficientEnumMapper[
                modelParameter.scalingCoefficientTypeEnum
            ]
            nonDimVal = func(modelParameter.dimensionalValue)
            modelParameter.setNondimensionalValue(nonDimVal)

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
