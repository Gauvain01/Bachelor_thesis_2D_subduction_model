import os

import numpy as np
from underworld import function as fn
from underworld import underworld as uw
from underworld import visualisation as vis

xRes = 192
yRes = 48
boxLength = 4.0
boxHeight = 1.0

mesh = uw.mesh.FeMesh_Cartesian(
    elementType=("Q1/DQ0"),
    elementRes=(xRes, yRes),
    maxCoord=(boxLength, boxHeight),
    periodic=[True, False],
)

velocityField = mesh.add_variable(nodeDofCount=2)
pressureField = mesh.add_variable(nodeDofCount=1)

swarm = uw.swarm.Swarm(mesh)

materialVariable = swarm.add_variable(dataType="int", count=1)
swarmLayout = uw.swarm.layouts.PerCellSpaceFillerLayout(swarm, 20)
swarm.populate_using_layout(swarmLayout)


# initialise the 'materialVariable' data to represent two different materials.
upperMantleIndex = 0
lowerMantleIndex = 1
upperSlabIndex = 2
lowerSlabIndex = 3
coreSlabIndex = 4
# Initial material layout has a flat lying slab with at 15\degree perturbation
lowerMantleY = 0.4
slabLowerShape = np.array(
    [
        (1.2, 0.925),
        (3.25, 0.925),
        (3.20, 0.900),
        (1.2, 0.900),
        (1.02, 0.825),
        (1.02, 0.850),
    ]
)
slabCoreShape = np.array(
    [
        (1.2, 0.975),
        (3.35, 0.975),
        (3.25, 0.925),
        (1.2, 0.925),
        (1.02, 0.850),
        (1.02, 0.900),
    ]
)
slabUpperShape = np.array(
    [
        (1.2, 1.000),
        (3.40, 1.000),
        (3.35, 0.975),
        (1.2, 0.975),
        (1.02, 0.900),
        (1.02, 0.925),
    ]
)

slabLower = uw.function.shape.Polygon(slabLowerShape)
slabCore = uw.function.shape.Polygon(slabCoreShape)
slabUpper = uw.function.shape.Polygon(slabUpperShape)

# initialise everying to be upper mantle material
materialVariable.data[:] = upperMantleIndex

lowerMantleY = 0.4

# change matieral index if the particle is not upper mantle
for index in range(len(swarm.particleCoordinates.data)):
    coord = swarm.particleCoordinates.data[index][:]
    if coord[1] < lowerMantleY:
        materialVariable.data[index] = lowerMantleIndex
    if slabCore.evaluate(tuple(coord)):
        materialVariable.data[index] = coreSlabIndex
    if slabUpper.evaluate(tuple(coord)):
        materialVariable.data[index] = upperSlabIndex
    elif slabLower.evaluate(tuple(coord)):
        materialVariable.data[index] = lowerSlabIndex

fig = vis.Figure(figsize=(960, 300), name="particles")
fig.append(
    vis.objects.Points(
        swarm, materialVariable, pointSize=2, colours="white green red purple blue"
    )
)
fig.show()
