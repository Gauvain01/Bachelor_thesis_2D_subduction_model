from distutils import core
import math


from underworld import function as fn
from underworld import mesh, scaling


from model_parameters.Model_parameter_set import ModelParameterSet

u = scaling.units


class RheologyCalculations:
    def __init__(
        self, modelParameterSet: ModelParameterSet, velocityField: mesh.MeshVariable
    ) -> None:
        self.modelParameterSet = modelParameterSet
        self.velocityField = velocityField

    def getSymmetricStrainRateTensor(self):
        symStrainRate = fn.tensor.symmetric(self.velocityField.fn_gradient)
        return symStrainRate

    def getStrainRateSecondInvariant(self):
        strainRateSecondInvariant = fn.tensor.second_invariant(
            fn.tensor.symmetric(self.velocityField.fn_gradient)
        )
        return strainRateSecondInvariant

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
        dt_e = self.modelParameterSet._scalingCoefficient.scalingForTime(2e4 * u.years)
        effVis = (coreVis * dt_e) / (alpha + dt_e)
        return effVis
