import math
from distutils import core

from underworld import function as fn
from underworld import mesh, mpi, scaling

from modelParameters import ModelParameterMap

u = scaling.units


class RheologyFunctions:
    def __init__(self, modelParameterMap: ModelParameterMap) -> None:
        self.modelParameterMap = modelParameterMap

        self.symStrainRate = None
        self.strainRateSecondInvariant = None
        self.rayLeighNumber = None

    def getSymmetricStrainRateTensor(self, velocityField):
        symStrainRate = fn.tensor.symmetric(velocityField.fn_gradient)
        return symStrainRate

    def getStrainRateSecondInvariant(self, velocityField):
        strainRateSecondInvariant = fn.tensor.second_invariant(
            fn.tensor.symmetric(velocityField.fn_gradient)
        )
        return strainRateSecondInvariant

    def getEffectiveViscosityOfUpperLayerVonMises(self, velocityField):

        sigmaY = (
            self.modelParameterMap.yieldStressOfSpTopLayer.nonDimensionalValue.magnitude
        )
        strainRateSecondInvariant = self.getStrainRateSecondInvariant(velocityField)

        effectiveViscosity = 0.5 * (sigmaY / (strainRateSecondInvariant + +1.0e-18))

        return effectiveViscosity

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
