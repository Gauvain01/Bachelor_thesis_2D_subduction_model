from enum import IntEnum


class ScalingCoefficientTypeEnum(IntEnum):
    VISCOSITY = 1
    STRESS = 2
    LENGTH = 3
    TEMPERATURE = 4
    GRADIENT = 5
    VELOCITYVECTOR = 6
    NONE = 7
    MASS = 8
    TIME = 9
    MANUAL = 10
