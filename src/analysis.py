import os
import pprint

from ModelAnalysis import ModelData, TracerAnalysis, TracerData
from strakParam import get_Strak_2021_model_parameter_map

# if __name__ is "__main__":
params = get_Strak_2021_model_parameter_map()
obj = ModelData("./src/output/test_69")
obj._getData()
tracer = obj.tracerDataSequence[2]
item = TracerAnalysis(tracer, params)
item.plotLateralVelocity("./src/output/upper_slabVel_69")
