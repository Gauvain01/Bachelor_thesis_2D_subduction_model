import attr

from modelParameters._Model_parameter import ModelParameter
from modelParameters._Scaling_coefficient import ScalingCoefficient


@attr.s(frozen=True, repr=True, slots=True, kw_only=True)
class ModelParameterMap:
    scalingCoefficient: ScalingCoefficient = attr.ib(
        validator=attr.validators.instance_of(ScalingCoefficient), repr=False
    )
    temperatureContrast: ModelParameter = attr.ib(
        validator=attr.validators.instance_of(ModelParameter)
    )
    modelHeight: ModelParameter = attr.ib(
        validator=attr.validators.instance_of(ModelParameter)
    )
    modelLength: ModelParameter = attr.ib(
        validator=attr.validators.instance_of(ModelParameter)
    )
    referenceDensity: ModelParameter = attr.ib(
        validator=attr.validators.instance_of(ModelParameter)
    )
    gravitationalAcceleration: ModelParameter = attr.ib(
        validator=attr.validators.instance_of(ModelParameter)
    )
    referenceTemperature: ModelParameter = attr.ib(
        validator=attr.validators.instance_of(ModelParameter)
    )
    thermalExpansivity: ModelParameter = attr.ib(
        validator=attr.validators.instance_of(ModelParameter)
    )
    thermalDiffusivity: ModelParameter = attr.ib(
        validator=attr.validators.instance_of(ModelParameter)
    )
    referenceViscosity: ModelParameter = attr.ib(
        validator=attr.validators.instance_of(ModelParameter)
    )
    upperMantleViscosity: ModelParameter = attr.ib(
        validator=attr.validators.instance_of(ModelParameter)
    )
    lowerMantleViscosity: ModelParameter = attr.ib(
        validator=attr.validators.instance_of(ModelParameter)
    )
    spTopLayerViscosity: ModelParameter = attr.ib(
        validator=attr.validators.instance_of(ModelParameter)
    )
    spCoreLayerViscosity: ModelParameter = attr.ib(
        validator=attr.validators.instance_of(ModelParameter)
    )
    spBottomLayerViscosity: ModelParameter = attr.ib(
        validator=attr.validators.instance_of(ModelParameter)
    )
    yieldStressOfSpTopLayer: ModelParameter = attr.ib(
        validator=attr.validators.instance_of(ModelParameter)
    )
    gasConstant: ModelParameter = attr.ib(
        validator=attr.validators.instance_of(ModelParameter)
    )
    lowerMantleHeigth: ModelParameter = attr.ib(
        validator=attr.validators.instance_of(ModelParameter)
    )
    coreShearModulus: ModelParameter = attr.ib(
        validator=attr.validators.instance_of(ModelParameter)
    )
    timeScaleStress: ModelParameter = attr.ib(
        validator=attr.validators.instance_of(ModelParameter)
    )
