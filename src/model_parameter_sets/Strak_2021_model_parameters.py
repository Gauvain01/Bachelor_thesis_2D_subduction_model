import UWGeodynamics as geo
from model_parameters.Model_parameter import ModelParameter
from model_parameters.Model_parameter_set import ModelParameterSet
from model_parameters.Scaling_coefficient import ScalingCoefficient
from model_parameters.scaling_coefficient_type_enum import ScalingCoefficientTypeEnum

u = geo.UnitRegistry


class Strak2021ScalingCoefficient(ScalingCoefficient):
    @property
    def viscosityCoefficient(self):
        return 3.5e20 * u.pascal * u.second

    @property
    def lengthCoefficient(self):
        return 660e3 * u.meter

    @property
    def temperatureCoefficient(self):
        return 1573.15 * u.kelvin

    @property
    def timeCoefficient(self):
        lens2 = self.lengthCoefficient**2
        k = (1e-6 * u.meter**2) / u.second
        return lens2 / k

    @property
    def gradientCoefficient(self):
        coeff = 1 / self.lengthCoefficient
        return coeff

    @property
    def massCoefficient(self):
        return 1.0 * u.kilogram

    @property
    def velocityVectorCoefficient(self):
        k = (1e-6 * u.meter**2) / u.second
        coeff = k / self.lengthCoefficient
        return coeff

    @property
    def stressCoefficient(self):
        k = (1e-6 * u.meter**2) / u.second
        visK = self.viscosityCoefficient * k
        lens2 = self.lengthCoefficient**2
        return visK / lens2


strak2021ModelParameterList = [
    ModelParameter("modelHeight", 660e3 * u.meter, ScalingCoefficientTypeEnum.LENGTH),
    ModelParameter("modelLength", 11600e3 * u.meter, ScalingCoefficientTypeEnum.LENGTH),
    ModelParameter(
        "referenceDensity",
        3230 * u.kilogram / u.meter**3,
        ScalingCoefficientTypeEnum.MANUAL,
    ).setNondimensionalValue(u.Quantity(0.0)),
    ModelParameter(
        "spUmDensityContrast",
        -62985 * u.kilogram / u.meter**3,
        ScalingCoefficientTypeEnum.MANUAL,
    ).setNondimensionalValue(u.Quantity(1.0)),
    ModelParameter(
        "gravitanionalAcceleration",
        9.81 * u.meter / u.second,
        ScalingCoefficientTypeEnum.NONE,
    ),
    ModelParameter(
        "referenceTemperature",
        1573.15 * u.kelvin,
        ScalingCoefficientTypeEnum.TEMPERATURE,
    ),
    ModelParameter(
        "spWmTemperatureContrast",
        1573.15 * u.kelvin,
        ScalingCoefficientTypeEnum.TEMPERATURE,
    ),
    ModelParameter(
        "thermalExpansivity", 1e-5 / u.kelvin, ScalingCoefficientTypeEnum.NONE
    ),
    ModelParameter(
        "thermalDiffusivity",
        1e-6 * u.meter**2 / u.second,
        ScalingCoefficientTypeEnum.NONE,
    ),
    ModelParameter(
        "referenceViscosity",
        3.5e20 * u.pascal * u.second,
        ScalingCoefficientTypeEnum.VISCOSITY,
    ),
    ModelParameter(
        "lowerViscosityCutoff",
        3.5e19 * u.pascal * u.second,
        ScalingCoefficientTypeEnum.VISCOSITY,
    ),
    ModelParameter(
        "upperViscosityCutoff",
        3.5e20 * u.pascal * u.second,
        ScalingCoefficientTypeEnum.VISCOSITY,
    ),
    ModelParameter(
        "lowerMantleViscosity",
        3.5e22 * u.pascal * u.second,
        ScalingCoefficientTypeEnum.VISCOSITY,
    ),
    ModelParameter(
        "spTopLayerViscosity",
        3.5e23 * u.pascal * u.second,
        ScalingCoefficientTypeEnum.VISCOSITY,
    ),
    ModelParameter(
        "opCrustLayerViscosity",
        3.5e23 * u.pascal * u.second,
        ScalingCoefficientTypeEnum.VISCOSITY,
    ),
    ModelParameter(
        "spCoreLayerViscosity",
        3.5e23 * u.pascal * u.second,
        ScalingCoefficientTypeEnum.VISCOSITY,
    ),
    ModelParameter(
        "spBottomLayerViscosity",
        1.75e22 * u.pascal * u.second,
        ScalingCoefficientTypeEnum.VISCOSITY,
    ),
    ModelParameter(
        "spEclogitizedTopLayerViscosity",
        1.75e22 * u.pascal * u.second,
        ScalingCoefficientTypeEnum.VISCOSITY,
    ),
    ModelParameter(
        "opLithosphericMantleViscosityInForearcAndBackarc",
        1.4e23 * u.pascal * u.second,
        ScalingCoefficientTypeEnum.VISCOSITY,
    ),
    ModelParameter(
        "opLithosphericMantleViscosityInFarBackarc",
        7e23 * u.pascal * u.second,
        ScalingCoefficientTypeEnum.VISCOSITY,
    ),
    ModelParameter(
        "yieldStressOfSpTopLayer", 21e6 * u.pascal, ScalingCoefficientTypeEnum.STRESS
    ),
    ModelParameter("preExponentialFactor", 3e6, ScalingCoefficientTypeEnum.NONE),
    ModelParameter(
        "activationEnergyUpperMantle",
        530e3 * u.joule / (u.mol),
        ScalingCoefficientTypeEnum.NONE,
    ),
    ModelParameter(
        "gasConstant",
        8.3145 * u.joule / (u.mol * u.kelvin),
        ScalingCoefficientTypeEnum.NONE,
    ),
    ModelParameter(
        "activationEnergyInSlab", None, ScalingCoefficientTypeEnum.MANUAL
    ).setNondimensionalValue(u.Quantity(0.0)),
]

strak2021ModelParameterSet = ModelParameterSet(Strak2021ScalingCoefficient())
strak2021ModelParameterSet.addModelParameterFromList(strak2021ModelParameterList)
strak2021ModelParameterSet.nonDimensionalizeParameters()
