from test.resources.Strak_2021_model_params_resources import getNonDimensionals

import numpy as np
from model_parameter_sets.Strak_2021_model_parameters import strak2021ModelParameterSet
from model_parameters.Model_parameter import ModelParameter
from RheologyFunctions import RheologyFunctions

# def test_strak_2021_model_params_scaling():
#     dimensionals = getNonDimensionals()
#     parameterList = strak2021ModelParameterSet
#     parameterList.nonDimensionalizeParameters()
#     testingList = []
#     for param in strak2021ModelParameterSet._modelParametersList:
#         param: ModelParameter
#         if param.nonDimensionalValue == None:
#             testingList.append(None)
#         else:
#             testingList.append(round(param.nonDimensionalValue.magnitude, 1))
#     testingArray = np.array(testingList)

#     assert np.array_equal(dimensionals, testingArray, equal_nan=True)


def test_rayleigh_number():
    rheologyFn = RheologyFunctions(strak2021ModelParameterSet, None)
    rayleighNumber = rheologyFn.getRayleighNumber()
    print(rayleighNumber)
    assert format(rayleighNumber, ".1E") == format(3.5e7, ".1E")
