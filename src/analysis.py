import os
import pprint

from ModelAnalysis import DipDataAnalysis, ModelData, TracerAnalysis, TracerData
from strakParam import get_Strak_2021_model_parameter_map

# if __name__ is "__main__":
params = get_Strak_2021_model_parameter_map()
obj = ModelData("./src/output/test_92")
obj._getData()
# tracer = obj.tracerDataSequence[0]
# item = TracerAnalysis(tracer, params)
# item.plotLateralVelocity("./src/output/hinge_92")

dip = obj.dipDataSequece[0]
item = DipDataAnalysis(dip, params)
item.plotDip("./src/output/mid_dip_92")
