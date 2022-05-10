from underworld import mesh, mpi, swarm, utils


def meshProjector(
    meshVariable: mesh.MeshVariable, swarmField: swarm.SwarmVariable
) -> None:
    with mpi.call_pattern():
        solver = utils.MeshVariable_Projection(meshVariable, swarmField)
        solver.solve()
