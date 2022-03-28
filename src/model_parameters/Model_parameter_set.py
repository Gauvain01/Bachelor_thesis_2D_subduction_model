import inspect
from typing import List

from matplotlib.pyplot import get

from model_parameters.Model_parameter import ModelParameter
from model_parameters.Scaling_coefficient import ScalingCoefficient


class ModelParameterSet:
    def __init__(self, scalingCoefficient: ScalingCoefficient) -> None:
        self._scalingCoefficient = scalingCoefficient
        self._modelParametersList = []

    def nonDimensionalizeParameters(self):
        parameters1: List[ModelParameter] = [
            a for a in inspect.getmembers(self) if not callable(a[1])
        ]
        print(parameters1)
        parameters = [a for a in parameters1 if not a[0].startswith("_")]
        print(parameters)
        for parameter in parameters:
            self._scalingCoefficient.nonDimensionalizeModelParameter(parameter[1])

    def addModelParameter(self, modelParameter: ModelParameter) -> None:
        try:
            getattr(self, modelParameter.name)
            raise ValueError
        except AttributeError:
            setattr(self, modelParameter.name, modelParameter)
            self._modelParametersList.append(modelParameter)

    def addModelParameterFromList(
        self, ModelParameterlist: List[ModelParameter]
    ) -> None:
        for param in ModelParameterlist:
            self.addModelParameter(param)

    def _printModelParametersNonDimensionalValues(self):
        for param in self._modelParametersList:
            print(f"{param.name} -> {param.nonDimensionalValue}")
