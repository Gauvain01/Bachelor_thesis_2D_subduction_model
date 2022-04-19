import math
from distutils import core

from underworld import function as fn
from underworld import mesh, scaling

from modelParameters import ModelParameterMap

u = scaling.units


class RheologyFunctions:
    def __init__(
        self, modelParameterMap: ModelParameterMap, velocityField: mesh.MeshVariable
    ) -> None:
        self.modelParameterMap = modelParameterMap
        self.velocityField = velocityField
        self.symStrainRate = None
        self.strainRateSecondInvariant = None
        self.rayLeighNumber = None

    def getSymmetricStrainRateTensor(self):
        if self.symStrainRate is None:
            self.symStrainRate = fn.tensor.symmetric(self.velocityField.fn_gradient)
            return self.symStrainRate
        else:
            return self.symStrainRate

    def getStrainRateSecondInvariant(self):
        if self.strainRateSecondInvariant is None:
            self.strainRateSecondInvariant = fn.tensor.second_invariant(
                fn.tensor.symmetric(self.velocityField.fn_gradient)
            )
            return self.strainRateSecondInvariant
        else:
            return self.strainRateSecondInvariant

    def getEffectiveViscosityOfUpperLayerVonMises(self):
        # sigmaY = (
        #     self.modelParameterMap.yieldStressOfSpTopLayer.nonDimensionalValue.magnitude
        # )
        strainRateSecondInvariant = self.getStrainRateSecondInvariant()

        # effectiveViscosity = (sigmaY) / (2 * strainRateSecondInvariant)

        # return effectiveViscosity
        cohesion = 0.06
        vonMises = 0.5 * cohesion / (strainRateSecondInvariant + 1.0e-18)
        return vonMises

    def getEffectiveViscosityOfViscoElasticCore(self):
        coreShearModulus = (
            self.modelParameterMap.coreShearModulus.nonDimensionalValue.magnitude
        )
        coreVis = (
            self.modelParameterMap.spCoreLayerViscosity.nonDimensionalValue.magnitude
        )

        alpha = coreVis / coreShearModulus
        dt_e = self.modelParameterMap.timeScaleStress.nonDimensionalValue.magnitude
        effVis = (coreVis * dt_e) / (alpha + dt_e)
        return effVis

    def getRayleighNumber(self):
        if self.rayLeighNumber is None:
            ls = self.modelParameterMap.modelHeight.dimensionalValue.magnitude
            print(f"{ls=}")
            rhoRef = self.modelParameterMap.referenceDensity.dimensionalValue.magnitude

            print(f"{rhoRef=}")
            g = (
                self.modelParameterMap.gravitationalAcceleration.dimensionalValue.magnitude
            )

            print(f"{g=}")
            alpha = self.modelParameterMap.thermalExpansivity.dimensionalValue.magnitude

            print(f"{alpha=}")
            deltaT = (
                self.modelParameterMap.temperatureContrast.dimensionalValue.magnitude
            )
            print(f"{deltaT=}")

            k = self.modelParameterMap.thermalDiffusivity.dimensionalValue.magnitude
            print(f"{k=}")

            eta = self.modelParameterMap.referenceViscosity.dimensionalValue.magnitude
            print(f"{eta=}")

            self.rayLeighNumber = ((ls**3) * rhoRef * g * alpha * deltaT) / (k * eta)
            return self.rayLeighNumber
        else:
            return self.rayLeighNumber
