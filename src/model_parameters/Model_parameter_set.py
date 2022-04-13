import inspect
from typing import List

from matplotlib.pyplot import get

from model_parameters.Model_parameter import ModelParameter
from model_parameters.Scaling_coefficient import ScalingCoefficient


class ModelParameterSet:
    def __init__(self, scalingCoefficient: ScalingCoefficient) -> None:
        self._scalingCoefficient = scalingCoefficient
        self._modelParametersList = []

        self.modelHeight: ModelParameter
        self.modelLength: ModelParameter
        self.referenceDensity: ModelParameter
        self.gravitationalAcceleration: ModelParameter
        self.referenceTemperature: ModelParameter
        self.thermalExpansivity: ModelParameter
        self.thermalDiffusivity: ModelParameter
        self.referenceViscosity: ModelParameter
        self.upperMantleViscosity: ModelParameter
        self.lowerMantleViscosity: ModelParameter
        self.spTopLayerViscosity: ModelParameter
        self.spCoreLayerViscosity: ModelParameter
        self.spBottomLayerViscosity: ModelParameter
        self.yieldStressOfSpTopLayer: ModelParameter
        self.gasConstant: ModelParameter
        self.lowerMantleHeigth: ModelParameter
        self.coreShearModulus: ModelParameter
        self.timeScaleStress: ModelParameter
        self._checked = False

    def _check(self):
        if not self._checked:
            parameters1: List[ModelParameter] = [
                a for a in inspect.getmembers(self) if not callable(a[1])
            ]
            parameters = [a for a in parameters1 if not a[0].startswith("_")]
            for parameter in parameters:
                if parameter is None:
                    raise ValueError(
                        f"{parameter} is None, add modelParameter with the correct name"
                    )
            self._checked = True

    def nonDimensionalizeParameters(self):
        self._check()
        parameters1: List[ModelParameter] = [
            a for a in inspect.getmembers(self) if not callable(a[1])
        ]
        parameters = [a for a in parameters1 if not a[0].startswith("_")]
        for parameter in parameters:
            self._scalingCoefficient.nonDimensionalizeModelParameter(parameter[1])

    def addModelParameter(self, modelParameter: ModelParameter) -> None:
        try:
            item = getattr(self, modelParameter.name)
            if item is not None:
                raise ValueError
            else:
                setattr(self, modelParameter.name, modelParameter)
        except AttributeError:
            setattr(self, modelParameter.name, modelParameter)
            self._modelParametersList.append(modelParameter)

    def addModelParameterFromList(
        self, ModelParameterlist: List[ModelParameter]
    ) -> None:
        for param in ModelParameterlist:
            self.addModelParameter(param)
        self._check()

    def _printModelParametersNonDimensionalValues(self):
        for param in self._modelParametersList:
            print(f"{param.name} -> {param.nonDimensionalValue}")
