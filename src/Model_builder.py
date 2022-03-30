import numpy as np
import UWGeodynamics as geo

from model_parameter_sets.Strak_2021_model_parameters import strak2021ModelParameterSet

parameterSet = strak2021ModelParameterSet
# nonDimensionalizeParameters
parameterSet.nonDimensionalizeParameters()
import math

import underworld as uw
import underworld.visualisation as vis

mesh = uw.mesh.FeMesh_Cartesian(
    elementRes=(4, 2),
    minCoord=(0.0, 0.0),
    maxCoord=(
        parameterSet.modelLength.nonDimensionalValue.magnitude,
        parameterSet.modelHeight.nonDimensionalValue.magnitude,
    ),
)

# visualising the result
figMesh = vis.Figure(figsize=(1080, 720))
figMesh.append(vis.objects.Mesh(mesh, nodeNumbers=True, pointsize=10))
figMesh.show()

# # setupModel
# Model = geo.Model(
#     elementRes=(1024, 512),
#     minCoord=(0, 0),
#     maxCoord=(
#         parameterSet.modelLength.dimensionValue,
#         parameterSet.modelHeight.dimensionValue,
#     ),
#     name="Test_1_With_Strak_2021_params",
#     outputDir="./output",
#     periodic=[False, False],
# )

# # let's define the shapes of the lower mantle
# upperMantleShape = geo.shapes.Layer(Model.top, Model.bottom)
# upperSlabShape = geo.shapes.Polygon(
#     np.array((530, 100), (5000, 563.04), (), (), (), ())
# )
