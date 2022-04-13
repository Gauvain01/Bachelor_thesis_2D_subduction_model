import math
from distutils import core

from underworld import function as fn
from underworld import mesh, scaling

from model_parameters.Model_parameter_set import ModelParameterSet

u = scaling.units


class RheologyFunctions:
    def __init__(
        self, modelParameterSet: ModelParameterSet, velocityField: mesh.MeshVariable
    ) -> None:
        self.modelParameterSet = modelParameterSet
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

    def getEffectiveViscosityDislocationCreep(self):
        # A * exp(E/nRT) * e_dot ^1-n/n
        E = (
            self.modelParameterSet.activationEnergyDislocationCreep.dimensionalValue.magnitude
        )
        A = (
            self.modelParameterSet.preExponentialFactorDislocation.dimensionalValue.magnitude
        )
        n = self.modelParameterSet.powerLawExponenet.nonDimensionalValue.magnitude
        R = self.modelParameterSet.gasConstant.dimensionalValue.magnitude
        T = self.modelParameterSet.referenceTemperature.dimensionalValue.magnitude
        strainRateSecondInvariant = self.getStrainRateSecondInvariant()

        dislocationViscosity = (
            A * math.exp((E / (n * R * T))) * strainRateSecondInvariant ** ((1 - n) / n)
        )
        scaledDislocationViscosity = (
            self.modelParameterSet._scalingCoefficient.scalingForViscosity(
                dislocationViscosity
            )
        )
        return scaledDislocationViscosity

    def getEffectiveViscosityDiffusionCreep(self):
        # A * exp(Ediff/RT)

        A = (
            self.modelParameterSet.preExponentialFactorDiffusion.dimensionalValue.magnitude
        )
        E = (
            self.modelParameterSet.activationEnergyDiffusionCreep.dimensionalValue.magnitude
        )
        R = self.modelParameterSet.gasConstant.dimensionalValue.magnitude
        T = self.modelParameterSet.referenceTemperature.dimensionalValue.magnitude

        diffusionViscosity = A * math.exp((E / (R * T)))
        scaledDiffusionViscosity = (
            self.modelParameterSet._scalingCoefficient.scalingForViscosity(
                diffusionViscosity
            )
        )
        return scaledDiffusionViscosity

    def getEffectiveViscosityOfUpperLayerVonMises(self):
        sigmaY = (
            self.modelParameterSet.yieldStressOfSpTopLayer.nonDimensionalValue.magnitude
        )
        strainRateSecondInvariant = self.getStrainRateSecondInvariant()

        effectiveViscosity = (sigmaY) / (2 * strainRateSecondInvariant)

        return effectiveViscosity

    def getEffectiveViscosityOfViscoElasticCore(self):
        coreShearModulus = self.modelParameterSet.coreShearModulus.nonDimensionalValue
        coreVis = self.modelParameterSet.spCoreLayerViscosity.nonDimensionalValue

        alpha = coreVis / coreShearModulus
        dt_e = self.modelParameterSet.timeScaleStress.nonDimensionalValue
        effVis = (coreVis * dt_e) / (alpha + dt_e)
        return effVis

    def getRayleighNumber(self):
        if self.rayLeighNumber is None:
            ls = self.modelParameterSet.modelHeight.dimensionalValue.magnitude
            print(ls)
            rhoRef = self.modelParameterSet.referenceDensity.dimensionalValue.magnitude
            print(rhoRef)
            g = (
                self.modelParameterSet.gravitationalAcceleration.dimensionalValue.magnitude
            )
            print(g)
            alpha = self.modelParameterSet.thermalExpansivity.dimensionalValue.magnitude
            print(alpha)
            deltaT = (
                self.modelParameterSet.temperatureContrast.dimensionalValue.magnitude
            )
            print(deltaT)
            k = self.modelParameterSet.thermalDiffusivity.dimensionalValue.magnitude
            print(k)
            eta = self.modelParameterSet.referenceViscosity.dimensionalValue.magnitude
            print(eta)

            self.rayLeighNumber = ((ls**3) * rhoRef * g * alpha * deltaT) / (k * eta)
            return self.rayLeighNumber
        else:
            return self.rayLeighNumber
