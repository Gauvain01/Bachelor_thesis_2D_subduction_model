from underworld import swarm, mesh, utils


def meshProjector(
    meshVariable: mesh.MeshVariable, swarmField: swarm.SwarmVariable
) -> None:
    solver = utils.MeshVariable_Projection(meshVariable, swarmField)
    solver.solve()
