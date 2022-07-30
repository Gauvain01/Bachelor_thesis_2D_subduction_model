import json

from underworld import function as fn
from underworld import mpi

from modelParameters import ModelParameterMap


class DipDataCollector:
    def __init__(
        self,
        minimumDepth: float,
        maximumDepth: float,
        modelParameters: ModelParameterMap,
        name: str,
    ) -> None:
        self.minimumDepth = minimumDepth
        self.maximumDepth = maximumDepth
        self.sca = modelParameters.scalingCoefficient
        self.params = modelParameters
        self.name = name
        if mpi.rank == 0:
            self.dipCoords = {}

    def setOutputPath(self, outputPath):
        self.outputPath = outputPath + f"/dip_{self.name}_coords.json"

    def calculateDip(
        self, materialVariable, swarm, materialIndex, modelTime, modelStep
    ):
        minimalDepth = (
            self.params.modelHeight.dimensionalValue.magnitude - self.minimumDepth
        ) / self.params.scalingCoefficient.lengthCoefficient.magnitude
        maximumDepth = (
            self.params.modelHeight.dimensionalValue.magnitude - self.maximumDepth
        ) / self.params.scalingCoefficient.lengthCoefficient.magnitude

        mapper = fn.branching.map(
            fn_key=materialVariable,
            mapping={materialIndex: fn.input()},
            fn_default=[0.0, 0.0],
        )
        item = mapper.evaluate(swarm)
        if mpi.rank == 0:
            totalCoordinateList = []
        with mpi.call_pattern("sequential"):
            coordinateList = []
            for i in item:
                if i[0] != 0.0 and i[1] != 0.0:
                    if i[1] >= minimalDepth and i[1] > maximumDepth:
                        if mpi.rank == 0:
                            totalCoordinateList.append([i[0], i[1]])
                        else:

                            coordinateList.append([i[0], i[1]])
            if mpi.rank != 0:
                mpi.comm.send(coordinateList, 0, tag=mpi.rank)

        mpi.barrier()
        if mpi.rank == 0:
            for i in range(1, mpi.size):
                coords = mpi.comm.recv(source=i, tag=i)
                if len(coords) != 0:
                    totalCoordinateList.append(coords)
                if len(totalCoordinateList) != 0:
                    coordDict = {"coord": totalCoordinateList, "time": modelTime}
                    self.dipCoords[modelStep] = coordDict

        mpi.barrier()

    def dipCoordsToJson(self):
        if mpi.rank == 0:
            path = f"{self.outputPath}/dip_shallow_coords_data.json"
            with open(path, "w") as f:
                json.dump(self.dipCoords, f)
        mpi.barrier()
