import math
from typing import List, Tuple

import numpy as np
import UWGeodynamics as geo
from UWGeodynamics import UnitRegistry

from model_parameters.Model_parameter_set import ModelParameterSet


class PlatePolygons:
    def __init__(
        self,
        parameterSet: ModelParameterSet,
        resolution: Tuple,
        dip: float,
        dipLength: UnitRegistry,
        plateLength: UnitRegistry,
        upperPlateThickness: UnitRegistry,
        middlePlateThickness: UnitRegistry,
        lowerPlateThickness: UnitRegistry,
        beginning: UnitRegistry,
    ) -> None:
        self.parameterSet = parameterSet
        self.parameterSet.nonDimensionalizeParameters()
        self.resolution = resolution
        self.dip = math.radians(dip)
        self._angle2 = 180 - 90 - dip
        self.otherDip = math.radians(self._angle2)
        self.dipLength = dipLength
        self.plateLength = plateLength
        self.u = geo.UnitRegistry
        self.upperPlateThickness = upperPlateThickness
        self.middlePlateThickness = middlePlateThickness
        self.lowerPlateThickness = lowerPlateThickness
        self.totalPlateThickness = (
            self.upperPlateThickness
            + self.middlePlateThickness
            + self.lowerPlateThickness
        )
        self.beginning = beginning
        self._calculatePolygons()

    def _calculatePolygons(self) -> np.ndarray:
        beginning = self.parameterSet._scalingCoefficient.scalingForLength(
            self.beginning
        ).magnitude
        print(beginning)

        lb = self.parameterSet._scalingCoefficient.scalingForLength(
            self.plateLength + self.beginning
        ).magnitude
        print(lb)
        dipLength = self.parameterSet._scalingCoefficient.scalingForLength(
            self.dipLength
        ).magnitude
        l1 = dipLength * math.cos(self.dip)
        print(l1)
        h1 = dipLength * math.sin(self.dip)
        print(h1)
        modelHeight = self.parameterSet.modelHeight.nonDimensionalValue.magnitude
        lowerThickness = self.parameterSet._scalingCoefficient.scalingForLength(
            self.lowerPlateThickness
        ).magnitude
        middleThickness = self.parameterSet._scalingCoefficient.scalingForLength(
            self.middlePlateThickness
        ).magnitude
        upperThickness = self.parameterSet._scalingCoefficient.scalingForLength(
            self.upperPlateThickness
        ).magnitude
        totalThick = lowerThickness + upperThickness + middleThickness
        coord1 = (lb + l1, modelHeight - h1 - totalThick)
        coord2 = (
            coord1[0] + lowerThickness * math.cos(self.otherDip),
            coord1[1] + lowerThickness * math.sin(self.otherDip),
        )
        coord3 = (
            coord2[0] + middleThickness * math.cos(self.otherDip),
            coord2[1] + middleThickness * math.sin(self.otherDip),
        )
        coord4 = (
            coord3[0] + upperThickness * math.cos(self.otherDip),
            coord3[1] + upperThickness * math.sin(self.otherDip),
        )

        coord8 = (lb, modelHeight - totalThick)
        coord7 = (
            coord8[0] + lowerThickness * math.cos(self.otherDip),
            coord8[1] + lowerThickness * math.sin(self.otherDip),
        )
        coord6 = (
            coord7[0] + middleThickness * math.cos(self.otherDip),
            coord7[1] + middleThickness * math.sin(self.otherDip),
        )
        coord5 = (
            coord6[0] + upperThickness * math.cos(self.otherDip),
            coord6[1] + upperThickness * math.sin(self.otherDip),
        )

        coord9 = (lb, modelHeight)
        coord10 = (lb, coord9[1] - upperThickness)
        coord11 = (lb, coord10[1] - middleThickness)

        coord13 = (beginning, modelHeight)
        coord14 = (beginning, coord10[1])
        coord15 = (beginning, coord11[1])
        coord16 = (beginning, coord8[1])

        self.lowerSlabPolygon = [
            coord15,
            coord11,
            coord7,
            coord2,
            coord1,
            coord8,
            coord16,
        ]

        self.middleSlabPolygon = [
            coord14,
            coord10,
            coord6,
            coord3,
            coord2,
            coord7,
            coord11,
            coord15,
        ]

        self.upperSlabPolygon = [
            coord13,
            coord9,
            coord5,
            coord4,
            coord3,
            coord6,
            coord10,
            coord14,
        ]

    def getUpperSlabShapeArray(self) -> List[Tuple]:
        return self.upperSlabPolygon

    def getMiddleSlabShapeArray(self) -> List[Tuple]:
        return self.middleSlabPolygon

    def getLowerSlabShapeArray(self) -> List[Tuple]:
        return self.lowerSlabPolygon
