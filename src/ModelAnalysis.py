import json
import os
from cProfile import label
from dataclasses import dataclass
from tracemalloc import start
from typing import Dict, List

import numpy as np
from matplotlib import pyplot as plt

from modelParameters import ModelParameterMap


class modelDataFolder:
    def __init__(self, path) -> None:
        self.path = path


class StepData(modelDataFolder):
    def __init__(self, path, stepNumber) -> None:
        super().__init__(path)
        self.h5Path = f"{self.path}/h5"
        self.xdmfPath = f"{self.path}/xdmf"
        self.timePath = f"{self.path}/time.json"
        stepNumber = stepNumber

    def getTimeJson(self) -> Dict:
        with open(self.timePath, "r") as f:
            time = json.loads(f)
        return time


class TracerData(modelDataFolder):
    def __init__(self, path) -> None:
        super().__init__(path)

    def getTracerData(self) -> Dict:
        with open(self.path, "r") as f:
            item = json.loads(f)
        return item


class ModelData:
    def __init__(self, dataPath) -> None:
        self.dataPath = dataPath
        self.stepDataSequence: List[StepData] = []
        self.tracerDataSequence: List[TracerData] = []

    def _getData(self):
        files = os.listdir(self.dataPath)

        for item in files:
            itemPath = f"{self.dataPath}/item"
            data = 1
            if len(item) == 5:
                try:
                    step = int(item)
                    folder = StepData(itemPath, step)
                    self.stepDataSequence.append(folder)
                    data = 0
                except:
                    ValueError
            if item.startswith("tracer"):
                tracerData = TracerData(itemPath)
                self.tracerDataSequence.append(tracerData)
                data = 0
            if item.startswith("mesh"):
                self.meshPath = itemPath

            if data == 1:
                raise ValueError(f"unaccounted item {item}")


class TracerAnalysis:
    def __init__(
        self, tracerData: TracerData, modelParameterMap: ModelParameterMap
    ) -> None:
        self.tracerData = tracerData
        self.modelParameters = modelParameterMap

    def plotLateralVelocity(self):
        lateralV, verticalV = self._calculateAverageVelocity()
        fig, ax = plt.subplots(figsize=(5, 2.7), layout="constrained")
        ax.plot([i[0] for i in lateralV], [i[1] for i in lateralV])
        ax.set_xlabel("Velocity (cm/yr)")
        ax.set_ylabel("Time (Myr)")
        ax.set_title("Lateral Velocity")
        fig.savefig("/src/output")

    def _calculateAverageVelocity(self):
        data = self.tracerData.getTracerData()
        startX = None
        startY = None
        lateralVelocity = []
        verticalVelocity = []
        for i in len(data):
            scaL = self.modelParameters.scalingCoefficient.lengthCoefficient.magnitude
            scaT = self.modelParameters.scalingCoefficient.timeCoefficient.magnitude
            stepData = data[str(i)]
            x = stepData["coord"][0] * scaL * 100
            y = stepData["coord"][1] * scaL * 100
            time = stepData["time"] * scaT / 60 / 60 / 24 / 365 / 1e6

            if i == 0:
                startX = x
                startY = y
                lateralVelocity.append((0.0, 0.0))
            if i > 0:
                deltaT = time - (
                    data[str(i) - 1]["time"] * scaT / 60 / 60 / 24 / 365 / 1e6
                )
                deltaX = abs(x - (data[str(i) - 1]["coord"][0] * scaL * 100))
                deltaY = abs(y - (data[str(i) - 1]["coord"][1] * scaL * 100))
                vVer = deltaY / deltaT
                vLat = deltaX / deltaT
                if y < startY:
                    vVer = -vVer
                if x < startX:
                    vlat = -vlat

                lateralVelocity.append((time, vLat))
                verticalVelocity.append((time, vVer))

        return lateralVelocity, verticalVelocity
