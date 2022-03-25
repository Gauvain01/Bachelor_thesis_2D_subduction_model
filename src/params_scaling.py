from dataclasses import dataclass

import numpy as np
from UWGeodynamics import UnitRegistry as u


@dataclass(frozen=True)
class MainPhysicalParameters:
    refDensity = 3300.0 * u.kilogram / u.meter**3
    refGravity = 9.8 * u.meter / u.second**2
    refDiffusivity = 1e-6 * u.meter**2 / u.second
    refExpansivity = 3e-5 / u.kelvin
    refViscosity = 1e20 * u.pascal * u.second
    refLength = 2900 * u.kilometer
    gasConstant = 8.314 * u.joule/(u.kilogram*u.kelvin)
    specificHeat = 1250.4 * u.joule/(u.kilogram * u.kelvin)
    potentialMantleTemp = 1673.*u.kelvin
    potentialTemp = 1673. * u.kelvin
    surfaceTemp = 273.*u.kelvin

@dataclass(frozen=True)
class MainRheologyParameters:
    cohesionMantle = 20.*u.megapascal
    frictionMantle = u.Quantity(0.1)
    frictionMantleDepth = frictionMantle*MainPhysicalParameters.refDensity*MainPhysicalParameters.refGravity
    diffusionPreExp = 1.87e9 * u.pascal*u.second
    diffusionEnergy = 3.16e5 * u.joule/(u.mol)
    diffusionEnergyDepth = diffusionEnergy * (1./MainPhysicalParameters.gasConstant)
    diffusionEnergyLowerMantle = 2e5 * u.joule/(u.mol)
    diffusionEnergyLowerMantleDepth = diffusionEnergyLowerMantle *(1/MainPhysicalParameters.gasConstant)
    diffusionVolume = 5.27e-6 * u.meter**3/(u.mol)
    diffusionVolumeDepth = diffusionVolume.magnitude * MainPhysicalParameters.refDensity.magnitude * MainPhysicalParameters.refGravity.magnitude * \
        u.joule/(u.mol * MainPhysicalParameters.gasConstant * u.meter)
    diffusionVolumeLowerMantle = 1.5e-6 * u.meter**3/(u.mol)
    diffusionVolumeLowerMantleDepth = diffusionVolumeLowerMantle.magnitude * MainPhysicalParameters.refDensity.magnitude * MainPhysicalParameters.refGravity.magnitude * \
        u.joule/(u.mol * MainPhysicalParameters.gasConstant * u.meter)
    adiabaticTempGradient = (MainPhysicalParameters.refExpansivity * MainPhysicalParameters.refGravity * MainPhysicalParameters.potentialTemp)/MainPhysicalParameters.specificHeat
    yieldStressMax = 200*u.megapascal
    lowerMantleViscosityFactor = u.Quantity(5.0)
    
    def __init__(self) -> None:
        dE = self.diffusionEnergy.magnitude - self.diffusionEnergyLowerMantle.magnitude
        dV = self.diffusionVolume.magnitude - self.diffusionEnergyLowerMantle.magnitude
        adT = self.adiabaticTempGradient.magnitude*660e3
        denom = MainPhysicalParameters.gasConstant.to_base_units().magnitude*(MainPhysicalParameters.potentialTemp.magnitude + adT)
        fac = np.exp((dE + MainPhysicalParameters.refDensity.magnitude*MainPhysicalParameters.refGravity.magnitude*660e3*dV)/(denom))
        self.diffusionPreExpLM = self.diffusionPreExp*fac
    
    
        
    
    


