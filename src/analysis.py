import os

from ModelAnalysis import ModelData, TracerAnalysis, TracerData
from strakParam import get_Strak_2021_model_parameter_map

# if __name__ is "__main__":
params = get_Strak_2021_model_parameter_map()

obj = ModelData("./src/output/test_67")
obj._getData()
tracer = obj.tracerDataSequence[0]
item = TracerAnalysis(tracer, params)
item.plotLateralVelocity()
