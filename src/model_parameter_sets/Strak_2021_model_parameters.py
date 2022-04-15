import UWGeodynamics as geo
from modelParameters import (
    ModelParameterBuilder,
    ScalingCoefficient,
    ModelParameterDaoBuilder,
    ScalingCoefficientType,
)

u = geo.UnitRegistry


class Strak2021ScalingCoefficient(ScalingCoefficient):
    def __init__(self) -> None:
        super().__init__()

    @property
    def viscosityCoefficient(self):
        return 3.5e20 * u.pascal * u.second

    @property
    def lengthCoefficient(self):
        return 2900e3 * u.meter

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
    ModelParameter("modelHeight", 2900e3 * u.meter, ScalingCoefficientType.LENGTH),
    ModelParameter("modelLength", 11600e3 * u.meter, ScalingCoefficientType.LENGTH),
    ModelParameter(
        "referenceDensity",
        3230 * u.kilogram / u.meter**3,
        ScalingCoefficientType.MANUAL,
    ).setNondimensionalValue(u.Quantity(0.0)),
    ModelParameter(
        "spUmDensityContrast",
        -62985 * u.kilogram / u.meter**3,
        ScalingCoefficientType.MANUAL,
    ).setNondimensionalValue(u.Quantity(1.0)),
    ModelParameter(
        "gravitationalAcceleration",
        9.81 * u.meter / u.second,
        ScalingCoefficientType.NONE,
    ),
    ModelParameter(
        "referenceTemperature",
        1573.15 * u.kelvin,
        ScalingCoefficientType.TEMPERATURE,
    ),
    ModelParameter(
        "temperatureContrast",
        1573.15 * u.kelvin,
        ScalingCoefficientType.TEMPERATURE,
    ),
    ModelParameter("thermalExpansivity", 1e-5 / u.kelvin, ScalingCoefficientType.NONE),
    ModelParameter(
        "thermalDiffusivity",
        1e-6 * u.meter**2 / u.second,
        ScalingCoefficientType.NONE,
    ),
    ModelParameter(
        "referenceViscosity",
        3.5e20 * u.pascal * u.second,
        ScalingCoefficientType.VISCOSITY,
    ),
    ModelParameter(
        "lowerViscosityCutoff",
        3.5e19 * u.pascal * u.second,
        ScalingCoefficientType.VISCOSITY,
    ),
    ModelParameter(
        "upperViscosityCutoff",
        3.5e20 * u.pascal * u.second,
        ScalingCoefficientType.VISCOSITY,
    ),
    ModelParameter(
        "lowerMantleViscosity",
        3.5e22 * u.pascal * u.second,
        ScalingCoefficientType.VISCOSITY,
    ),
    ModelParameter(
        "spTopLayerViscosity",
        3.5e23 * u.pascal * u.second,
        ScalingCoefficientType.VISCOSITY,
    ),
    ModelParameter(
        "opCrustLayerViscosity",
        3.5e23 * u.pascal * u.second,
        ScalingCoefficientType.VISCOSITY,
    ),
    ModelParameter(
        "spCoreLayerViscosity",
        3.5e23 * u.pascal * u.second,
        ScalingCoefficientType.VISCOSITY,
    ),
    ModelParameter(
        "spBottomLayerViscosity",
        1.75e22 * u.pascal * u.second,
        ScalingCoefficientType.VISCOSITY,
    ),
    ModelParameter(
        "spEclogitizedTopLayerViscosity",
        1.75e22 * u.pascal * u.second,
        ScalingCoefficientType.VISCOSITY,
    ),
    ModelParameter(
        "opLithosphericMantleViscosityInForearcAndBackarc",
        1.4e23 * u.pascal * u.second,
        ScalingCoefficientType.VISCOSITY,
    ),
    ModelParameter(
        "opLithosphericMantleViscosityInFarBackarc",
        7e23 * u.pascal * u.second,
        ScalingCoefficientType.VISCOSITY,
    ),
    ModelParameter(
        "yieldStressOfSpTopLayer", 21e6 * u.pascal, ScalingCoefficientType.STRESS
    ),
    ModelParameter("preExponentialFactor", 3e6, ScalingCoefficientType.NONE),
    ModelParameter(
        "activationEnergyUpperMantle",
        530e3 * u.joule / (u.mol),
        ScalingCoefficientType.NONE,
    ),
    ModelParameter(
        "gasConstant",
        8.3145 * u.joule / (u.mol * u.kelvin),
        ScalingCoefficientType.NONE,
    ),
    ModelParameter(
        "activationEnergyInSlab", None, ScalingCoefficientType.MANUAL
    ).setNondimensionalValue(u.Quantity(0.0)),
]

modelParameterBuilder = ModelParameterBuilder(Strak2021ScalingCoefficient())
builder = ModelParameterDaoBuilder(modelParameterBuilder)

strakParameterDao = (
    builder.setGasConstant(u.Quantity(8.3145))
    .setYieldStressOfSpTopLayer(21e6 * u.pascal, ScalingCoefficientType.STRESS)
    .setCoreShearModulus(u.Quantity(1e4))
    .setLowerMantleHeight(660e3 * u.meter, ScalingCoefficientType.LENGTH)
    .setModelHeight(2900e3 * u.meter, ScalingCoefficientType.LENGTH)
    .setModelLength(11600e3 * u.meter, ScalingCoefficientType.LENGTH)
    .setGravitationalAcceleration(
        9.81 * u.meter / u.second**2, ScalingCoefficientType.NONE
    )
    .setLowerMantleViscosity(
        3.5e22 * u.pascal * u.second, ScalingCoefficientType.VISCOSITY
    )
    .setUpperMantleViscosity(
        3.5e20 * u.pascal * u.second, ScalingCoefficientType.VISCOSITY
    )
    .set
)
