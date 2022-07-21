import UWGeodynamics as geo
from modelParameters import (
    ModelParameterBuilder,
    ModelParameterDaoBuilder,
    ScalingCoefficient,
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


modelParameterBuilder = ModelParameterBuilder(Strak2021ScalingCoefficient())
builder = ModelParameterDaoBuilder(modelParameterBuilder)

StrakParameterDao = (
    builder.setGasConstant(u.Quantity(8.3145))
    .setYieldStressOfSpTopLayer(21e6 * u.pascal, ScalingCoefficientType.STRESS)
    .setCoreShearModulus(u.Quantity(1e4))
    .setLowerMantleHeight(660e3 * u.meter, ScalingCoefficientType.LENGTH)
    .setModelHeight(660e3 * u.meter, ScalingCoefficientType.LENGTH)
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
    .setReferenceDensity(u.Quantity(0.0))
    .setThermalDiffusivity(1e-6 * u.meter**2 / u.second, ScalingCoefficientType.NONE)
    .setThermalExpansivity(1e-5 * u.meter**2 / u.second, ScalingCoefficientType.NONE)
    .setReferenceTemperature(1573.15 * u.kelvin, ScalingCoefficientType.TEMPERATURE)
    .setReferenceViscosity(
        3.5e20 * u.pascal * u.second, ScalingCoefficientType.VISCOSITY
    )
    .setSpBottomLayerViscosity(
        1.75e22 * u.pascal * u.second, ScalingCoefficientType.VISCOSITY
    )
    .setSpCoreLayerViscosity(
        3.5e23 * u.pascal * u.second, ScalingCoefficientType.VISCOSITY
    )
    .setSpTopLayerViscosity(
        3.5e23 * u.pascal * u.second, ScalingCoefficientType.VISCOSITY
    )
    .setTimeScaleStress(u.Quantity(1.0))
)
print(StrakParameterDao.__dict__)
