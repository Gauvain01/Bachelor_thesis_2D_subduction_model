import attr

from modelParameters._Model_parameter import ModelParameter


@attr.s(frozen=True, repr=True, slots=True)
class ModelParameterDao:
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

    def __dict__(self):
        return {key: item.__dict__() for key, item in super().__dict__()}
