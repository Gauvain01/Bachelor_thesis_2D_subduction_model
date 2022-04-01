from __future__ import annotations

import inspect
from dataclasses import dataclass

import numpy as np
import UWGeodynamics as geo

u = geo.u


class NonDimensionalizer:
    def __init__(self) -> None:
        self._setScalingCoefficients()

    def _setScalingCoefficients(self):
        KL = MainPhysicalParameters.refLength
        KT = MainPhysicalParameters.potentialTemp - MainPhysicalParameters.surfaceTemp
        Kt = KL**2 / MainPhysicalParameters.refDiffusivity
        KM = MainPhysicalParameters.refViscosity * KL * Kt
        Kvis = MainPhysicalParameters.refViscosity
        Ksig = (
            MainPhysicalParameters.refDiffusivity * MainPhysicalParameters.refViscosity
        ) / KL**2
        Kgrad = 1 / KL
        Kvel = MainPhysicalParameters.refViscosity / KL

        geo.scaling_coefficients["[length]"] = KL
        geo.scaling_coefficients["[viscosity]"] = Kvis
        geo.scaling_coefficients["[temperature]"] = KT

        # geo.scaling_coefficients["[mass]"] = KM
        geo.scaling_coefficients["[time]"] = Kt
        geo.scaling_coefficients["[velocity]"] = Kvel

        geo.scaling_coefficients["[stress]"] = Ksig
        geo.scaling_coefficients["[gradient]"] = Kgrad


@dataclass(frozen=True)
class MainPhysicalParameters:
    refDensity = 3300.0 * u.kilogram / u.meter**3
    refGravity = 9.8 * u.meter / u.second**2
    refDiffusivity = 1e-6 * u.meter**2 / u.second
    refExpansivity = 3e-5 / u.kelvin
    refViscosity = 1e20 * u.pascal * u.second
    refLength = 2900 * u.kilometer
    gasConstant = 8.314 * u.joule / (u.kilogram * u.kelvin)
    specificHeat = 1250.4 * u.joule / (u.kilogram * u.kelvin)
    potentialMantleTemp = 1673.0 * u.kelvin
    potentialTemp = 1673.0 * u.kelvin
    surfaceTemp = 273.0 * u.kelvin

    def getNonDimensionalized(self) -> MainPhysicalParametersNonDimensionalized:
        return MainPhysicalParametersNonDimensionalized()


@dataclass(frozen=True)
class MainPhysicalParametersNonDimensionalized(NonDimensionalizer):
    def __init__(self) -> None:
        super().__init__()
        dimensionalValues = self._getDimensionalValues()
        for i in dimensionalValues:
            object.__setattr__(self, i[0], geo.non_dimensionalise(i[1]))

    def _getDimensionalValues(self):
        return [
            a
            for a in inspect.getmembers(MainPhysicalParameters())
            if not a[0].startswith("_")
        ]


@dataclass(frozen=True)
class MainRheologyParameters:
    cohesionMantle = 20.0 * u.megapascal
    frictionMantle = u.Quantity(0.1)
    frictionMantleDepth = (
        frictionMantle
        * MainPhysicalParameters.refDensity
        * MainPhysicalParameters.refGravity
    )
    diffusionPreExp = 1.87e9 * u.pascal * u.second
    diffusionEnergy = 3.16e5 * u.joule / (u.mol)
    diffusionEnergyDepth = diffusionEnergy * (1.0 / MainPhysicalParameters.gasConstant)
    diffusionEnergyLowerMantle = 2e5 * u.joule / (u.mol)
    diffusionEnergyLowerMantleDepth = diffusionEnergyLowerMantle * (
        1 / MainPhysicalParameters.gasConstant
    )
    diffusionVolume = 5.27e-6 * u.meter**3 / (u.mol)
    diffusionVolumeDepth = (
        diffusionVolume.magnitude
        * MainPhysicalParameters.refDensity.magnitude
        * MainPhysicalParameters.refGravity.magnitude
        * u.joule
        / (u.mol * MainPhysicalParameters.gasConstant * u.meter)
    )
    diffusionVolumeLowerMantle = 1.5e-6 * u.meter**3 / (u.mol)
    diffusionVolumeLowerMantleDepth = (
        diffusionVolumeLowerMantle.magnitude
        * MainPhysicalParameters.refDensity.magnitude
        * MainPhysicalParameters.refGravity.magnitude
        * u.joule
        / (u.mol * MainPhysicalParameters.gasConstant * u.meter)
    )
    adiabaticTempGradient = (
        MainPhysicalParameters.refExpansivity
        * MainPhysicalParameters.refGravity
        * MainPhysicalParameters.potentialTemp
    ) / MainPhysicalParameters.specificHeat
    yieldStressMax = 200 * u.megapascal
    lowerMantleViscosityFactor = u.Quantity(5.0)

    def __init__(self) -> None:
        dE = self.diffusionEnergy.magnitude - self.diffusionEnergyLowerMantle.magnitude
        dV = self.diffusionVolume.magnitude - self.diffusionVolumeLowerMantle.magnitude
        adT = self.adiabaticTempGradient.magnitude * 660e3
        denom = MainPhysicalParameters.gasConstant.to_base_units().magnitude * (
            MainPhysicalParameters.potentialTemp.magnitude + adT
        )
        fac = np.exp(
            (
                dE
                + MainPhysicalParameters.refDensity.magnitude
                * MainPhysicalParameters.refGravity.magnitude
                * 660e3
                * dV
            )
            / (denom)
        )
        diffusionPreExpLM = self.diffusionPreExp * fac
        object.__setattr__(self, "diffusionPreExpLM", diffusionPreExpLM)

    def getNonDimensionalized(self) -> MainPhysicalParametersNonDimensionalized:
        return MainRheologyParametersNonDimensionalized()


@dataclass(frozen=True)
class MainRheologyParametersNonDimensionalized(NonDimensionalizer):
    def __init__(self) -> None:
        super().__init__()
        dimensionalValues = self._getDimensionalValues()
        for i in dimensionalValues:
            object.__setattr__(self, i[0], geo.non_dimensionalise(i[1]))

    def _getDimensionalValues(self):
        return [
            a
            for a in inspect.getmembers(MainRheologyParameters())
            if not a[0].startswith("_")
        ]
