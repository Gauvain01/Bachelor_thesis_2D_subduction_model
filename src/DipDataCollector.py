import json
import math
from typing import Dict, Tuple

import numpy as np
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
            self.dipAngles = {}

        self.coordinateList = {}

    def setOutputPath(self, outputPath):
        self.outputPath = outputPath + f"/dip_{self.name}_angles.json"

    def calculateDip(
        self, materialVariable, swarm, materialIndex, modelTime, modelStep
    ):
        if mpi.rank == 0:
            lowestCoordsOuput = []
            highestCoordsOutput = []
        mpi.barrier()
        minimalDepth = (
            self.params.modelHeight.dimensionalValue.magnitude - self.minimumDepth
        ) / self.params.scalingCoefficient.lengthCoefficient.magnitude
        maximumDepth = (
            self.params.modelHeight.dimensionalValue.magnitude - self.maximumDepth
        ) / self.params.scalingCoefficient.lengthCoefficient.magnitude
        variableL = 1e3 / self.params.scalingCoefficient.lengthCoefficient.magnitude
        mapper = fn.branching.map(
            fn_key=materialVariable,
            mapping={materialIndex: fn.input()},
            fn_default=[0.0, 0.0],
        )
        item = mapper.evaluate(swarm)

        with mpi.call_pattern("sequential"):
            lowRangeCoords = {}
            lowRangeCoordsList = []
            highRangeCoords = {}
            highRangeCoordsList = []
            # getting the relevant materialIndex cooords
            for i in item:
                if i[0] != 0.0 and i[1] != 0.0:
                    # finding the cords with the correct depth
                    if i[1] < maximumDepth and i[1] > (maximumDepth - variableL):
                        lowRangeCoordsList.append([i[0], i[1]])

                    if i[1] < minimalDepth and i[1] > (minimalDepth - variableL):
                        highRangeCoordsList.append([i[0], i[1]])
            # sort for the cords that have the closest y value to the target, keys are x values
            if len(highRangeCoordsList) != 0:
                for i in highRangeCoordsList:
                    try:
                        coord = highRangeCoords[i[0]]
                        if i[1] > coord[1]:
                            highRangeCoords[i[0]] = i
                    except KeyError:
                        highRangeCoords[i[0]] = i
                # finding the highest X value for the keys
                highestCoordKeys = [i for i in highRangeCoords.keys()]
                highestCoordKeys.sort()
                highestCoordKeyWithHighX = highestCoordKeys.pop()
                highestCoord = highRangeCoords[highestCoordKeyWithHighX]
            else:
                highestCoord = None

            if len(lowRangeCoordsList) != 0:
                for i in lowRangeCoordsList:
                    try:
                        coord = lowRangeCoords[i[0]]
                        if i[1] > coord[1]:
                            lowRangeCoords[i[0]] = i
                    except KeyError:
                        lowRangeCoords[i[0]] = i
                # finding the highest X value for the keys
                sortedLowestCoordKes = [i for i in lowRangeCoords.keys()]
                sortedLowestCoordKes.sort()
                lowestCoordKeyHighX = sortedLowestCoordKes.pop()
                lowestCoord = lowRangeCoords[lowestCoordKeyHighX]
            else:
                lowestCoord = None

            output = (highestCoord, lowestCoord)
            if mpi.rank != 0:
                mpi.comm.send(obj=output, dest=0, tag=23)
            else:
                if output[1] is not None:
                    lowestCoordsOuput.append(output[1])
                if output[0] is not None:
                    highestCoordsOutput.append(output[0])

        mpi.barrier()
        if mpi.rank == 0:
            iterList = [i for i in range(mpi.size)]
            iterList.pop(0)
            for i in iterList:
                item: Tuple = mpi.comm.recv(source=i, tag=23)
                if item[0] is not None:
                    highestCoordsOutput.append(item[0])
                if item[1] is not None:
                    lowestCoordsOuput.append(item[1])
            if len(highestCoordsOutput) != 0 or len(lowestCoordsOuput) != 0:
                highestCoordsOutput.sort()
                lowestCoordsOuput.sort()
                finalLowestCoord = lowestCoordsOuput.pop()
                finalHighestCoord = highestCoordsOutput.pop()

                deltaHScaled = finalHighestCoord[1] - finalLowestCoord[1]
                deltaLScaled = abs(finalHighestCoord[0] - finalLowestCoord[0])

                deltaHDim = deltaHScaled * self.sca.lengthCoefficient.magnitude
                deltaLDim = deltaLScaled * self.sca.lengthCoefficient.magnitude
                angleRad = math.atan((deltaHDim / deltaLDim))
                angle = math.degrees(angleRad)
                self.dipAngles[modelStep] = {"angle": angle, "time": modelTime}
        mpi.barrier()

    def dipAnglesToJson(self):
        mpi.barrier()
        if mpi.rank == 0:
            path = self.outputPath
            with open(path, "w") as f:
                json.dump(self.dipAngles, f)
        mpi.barrier()
