from underworld import underworld as uw
from underworld import visualisation

myMesh = uw.mesh.FeMesh_Cartesian(
    elementRes=(8, 8), minCoord=(0.0, 0.0), maxCoord=(2.0, 1.0)
)

temperatureField = myMesh.add_variable(nodeDofCount=1)

for index, coord in enumerate(myMesh.data):
    temperatureField.data[index] = coord[1]

fig = visualisation.Figure(figsize=(800, 400))
fig.append(
    visualisation.objects.Surface(myMesh, temperatureField, colours="blue white red")
)
fig.append(visualisation.objects.Mesh(myMesh))
fig.show()
