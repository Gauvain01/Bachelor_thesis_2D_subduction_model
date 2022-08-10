import json
import os
from cProfile import label
from dataclasses import dataclass
from tracemalloc import start
from typing import Dict, List

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axis import Axis

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
            time = json.load(f)
        return time


class DipData(modelDataFolder):
    def __init__(self, path) -> None:
        super().__init__(path)

    def getDipData(self) -> Dict:
        with open(self.path, "r") as f:
            item = json.load(f)
        return item


class TracerData(modelDataFolder):
    def __init__(self, path) -> None:
        super().__init__(path)

    def getTracerData(self) -> Dict:
        with open(self.path, "r") as f:
            item = json.load(f)
        return item


class ModelData:
    def __init__(self, dataPath) -> None:
        self.dataPath = dataPath
        self.stepDataSequence: List[StepData] = []
        self.tracerDataSequence: List[TracerData] = []
        self.dipDataSequece: List[DipData] = []

    def _getData(self):
        files = os.listdir(self.dataPath)

        for item in files:
            itemPath = f"{self.dataPath}/{item}"
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
            if item.startswith("dip"):
                dipData = DipData(itemPath)
                self.dipDataSequece.append(dipData)
                data = 0
            if item.startswith("mesh") or item.endswith("FigStore.gldb"):
                self.meshPath = itemPath
                data = 0

            if data == 1:
                raise ValueError(f"unaccounted item {item}")


class TracerAnalysis:
    def __init__(
        self, tracerData: TracerData, modelParameterMap: ModelParameterMap
    ) -> None:
        self.tracerData = tracerData
        self.modelParameters = modelParameterMap

    def plotLateralVelocity(self, path):
        lateralV, verticalV = self._calculateAverageVelocity()
        fig, ax = plt.subplots(figsize=(5, 2.7), layout="constrained")
        ax.plot([i[0] for i in lateralV], [i[1] for i in lateralV])
        ax.set_ylabel("Velocity (cm/yr)")
        ax.set_xlabel("Time (Myr)")
        ax.set_title("Lateral Velocity")
        self.tracerData.path
        fig.savefig(path)

    def _calculateAverageVelocity(self):
        data = self.tracerData.getTracerData()
        startX = None
        startY = None
        lateralVelocity = []
        verticalVelocity = []
        steps = [int(i) for i in data.keys()]
        steps.sort()

        for i in steps:
            scaL = 2900e3
            k = self.modelParameters.thermalDiffusivity.dimensionalValue.magnitude

            lSquared = 2900e3
            scaT = (lSquared * lSquared) / k

            stepData = data[f"{i}"]
            x = stepData["coord"][0] * scaL * 100
            y = stepData["coord"][1] * scaL * 100
            time = (stepData["time"] * scaT) / (60 * 60 * 24 * 365 * 1e6)

            if i == 0:
                startX = x
                startY = y
                lateralVelocity.append((0.0, 0.0))
            if i > 0:
                deltaT = time - (
                    data[f"{i-1}"]["time"] * scaT / (60 * 60 * 24 * 365 * 1e6)
                )
                deltaX = abs(x - (data[f"{i-1}"]["coord"][0] * scaL * 100))
                deltaY = abs(y - (data[f"{i-1}"]["coord"][1] * scaL * 100))
                vVer = deltaY / (deltaT * 1e6)
                vLat = deltaX / (deltaT * 1e6)
                if y < startY:
                    vVer = -vVer
                if x < startX:
                    vLat = -vLat

                lateralVelocity.append((time, vLat))
                verticalVelocity.append((time, vVer))

        return lateralVelocity, verticalVelocity


class DipDataAnalysis:
    def __init__(self, dipData: DipData, modelParameterMap: ModelParameterMap) -> None:
        self.dipData = dipData
        self.modelParameters = modelParameterMap

    def calculateDip(self):
        data = self.dipData.getDipData()
        steps = [int(i) for i in data.keys()]
        steps.sort()
        timeDegreePair = []
        for i in steps:
            k = self.modelParameters.thermalDiffusivity.dimensionalValue.magnitude

            lSquared = 2900e3
            scaT = (lSquared * lSquared) / k

            stepData = data[f"{i}"]

            time = (stepData["time"] * scaT) / (60 * 60 * 24 * 365 * 1e6)
            angle = stepData["angle"]

            timeDegreePair.append((time, angle))

        return timeDegreePair

    def plotDip(self, path):
        timeDegreePair = self.calculateDip()
        timeList = [i[0] for i in timeDegreePair]
        angles = [i[1] for i in timeDegreePair]
        fig, ax = plt.subplots(figsize=(5, 2.7), layout="constrained")
        ax.set_yticks([i for i in range(-20, 100, 10)])
        ax.scatter(timeList, angles, marker="p")
        ax.set_yticks([i for i in range(-20, 100, 10)])
        ax.set_ylabel("dip angle(degrees)")
        ax.set_xlabel("Time (Myr)")
        ax.set_title("dip")
        fig.savefig(path)
