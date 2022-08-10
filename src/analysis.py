import os
import pprint

from ModelAnalysis import DipDataAnalysis, ModelData, TracerAnalysis, TracerData
from strakParam import get_Strak_2021_model_parameter_map

# if __name__ is "__main__":
params = get_Strak_2021_model_parameter_map()
obj = ModelData("./src/output/test_84")
obj._getData()
# tracer = obj.tracerDataSequence[2]
# item = TracerAnalysis(tracer, params)
# item.plotLateralVelocity("./src/output/slab_84")

dip = obj.dipDataSequece[1]
item = DipDataAnalysis(dip, params)
item.plotDip("./src/output/shallow_dip_84")
