import math
from typing import List, Tuple

import numpy as np
import UWGeodynamics as geo
from UWGeodynamics import UnitRegistry

from model_parameters.Model_parameter_set import ModelParameterSet


class OceanicPLatePolygons:
    def __init__(
        self,
        parameterSet: ModelParameterSet,
        dip: float,
        dipLength: UnitRegistry,
        plateLength: UnitRegistry,
        upperPlateThickness: UnitRegistry,
        middlePlateThickness: UnitRegistry,
        lowerPlateThickness: UnitRegistry,
        beginning: UnitRegistry,
        foreArcAndBackArcLength,
        transitionLength,
        farBackArcLength,
        crustThickness,
        lithosphericForeArcThickness,
        lithospehricFarBackArcThickness,
    ) -> None:
        self.parameterSet = parameterSet
        self.parameterSet.nonDimensionalizeParameters()
        self.dip = math.radians(dip)
        self._angle2 = 180 - 90 - dip
        self.otherDip = math.radians(self._angle2)
        self.dipLength = dipLength
        self.plateLength = plateLength
        self.u = geo.UnitRegistry
        self.upperPlateThickness = upperPlateThickness
        self.middlePlateThickness = middlePlateThickness
        self.lowerPlateThickness = lowerPlateThickness
        self.beginning = beginning
        self.foreArcAndBackArcLength = foreArcAndBackArcLength
        self.transitionLength = transitionLength
        self.farBackArcLength = farBackArcLength
        self.crustThickness = crustThickness
        self.lithosphericForeArcThickness = lithosphericForeArcThickness
        self.lithosphericFarBackArcThickness = lithospehricFarBackArcThickness

        self._calculatePolygons()

    def _calculatePolygons(self) -> np.ndarray:
        beginning = self.parameterSet._scalingCoefficient.scalingForLength(
            self.beginning
        ).magnitude

        lb = self.parameterSet._scalingCoefficient.scalingForLength(
            self.plateLength + self.beginning
        ).magnitude

        dipLength = self.parameterSet._scalingCoefficient.scalingForLength(
            self.dipLength
        ).magnitude
        l1 = dipLength * math.cos(self.dip)

        h1 = dipLength * math.sin(self.dip)

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
        self.coord5 = (
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
            self.coord5,
            coord4,
            coord3,
            coord6,
            coord10,
            coord14,
        ]

        # calculate function of equal legged triangle
        # y=ax+b
        self._magnitude = (
            self.parameterSet._scalingCoefficient.lengthCoefficient.magnitude
        )
        self.elt_a = (
            coord9[1] * self._magnitude - self.coord5[1] * self._magnitude
        ) / (coord9[0] * self._magnitude - self.coord5[0] * self._magnitude)
        # y - coord9[1] = self.elt_a * (x - coord9[0])
        # y = self.elt_a x - coord9[0] * self.elt_a + coord9[1]
        self.elt_b = (
            -self.elt_a * coord9[0] * self._magnitude + coord9[1] * self._magnitude
        )

        # calculate function of upperPart of slab
        # y=ax+b

        self.ups_a = (
            self.coord5[1] * self._magnitude - coord4[1] * self._magnitude
        ) / (self.coord5[0] * self._magnitude - coord4[0] * self._magnitude)
        self.ups_b = (
            -self.ups_a * coord4[0] * self._magnitude + coord4[1] * self._magnitude
        )

        foreArcAndBackArcLength = (
            self.parameterSet._scalingCoefficient.scalingForLength(
                self.foreArcAndBackArcLength
            ).magnitude
        )
        transitionLength = self.parameterSet._scalingCoefficient.scalingForLength(
            self.transitionLength
        ).magnitude
        farBackarclength = self.parameterSet._scalingCoefficient.scalingForLength(
            self.farBackArcLength
        ).magnitude
        crustThickness = self.parameterSet._scalingCoefficient.scalingForLength(
            self.crustThickness
        ).magnitude
        lithosphericForArcThickness = (
            self.parameterSet._scalingCoefficient.scalingForLength(
                self.lithosphericForeArcThickness
            ).magnitude
        )
        lithosphericFarBackArcThicknes = (
            self.parameterSet._scalingCoefficient.scalingForLength(
                self.lithosphericFarBackArcThickness
            ).magnitude
        )
        totalLenLithosphericPLate = (
            foreArcAndBackArcLength + transitionLength + farBackarclength
        )
        # coordinatesCrustLayer
        coord21 = (coord9[0] + totalLenLithosphericPLate, modelHeight)
        coord22 = (coord21[0], coord21[1] - crustThickness)
        if self._isNodeOnEqualTriangle(coord22[1]):
            print("4")
            coord23 = (self._equalLeggedTriangle(coord22[1]), coord22[1])
            self.crustSlabShape = [coord9, coord21, coord22, coord23]
        else:
            print("5")
            coord23 = (self._upsCoord(coord22[1]), coord22[1])
            self.crustSlabShape = [coord9, coord21, coord22, coord23, self.coord5]

        # coordinates lithosphericMantleLayer
        coord28 = (coord22[0], coord22[1] - lithosphericFarBackArcThicknes)
        coord25 = (coord22[0] - farBackarclength, coord28[1])
        coord26 = (
            coord25[0] - transitionLength,
            coord22[1] - lithosphericForArcThickness,
        )
        coord24 = (coord25[0], coord22[1])

        if self._isNodeOnEqualTriangle(coord26[1]):
            print("1")
            coord27 = (self._equalLeggedTriangle(coord26[1]), coord26[1])
            self.lithosSphericMantleShape = [
                coord23,
                coord24,
                coord25,
                coord26,
                coord27,
            ]
        if self._isNodeOnEqualTriangle(
            coord26[1]
        ) == False and self._isNodeOnEqualTriangle(coord23[1]):
            print("2")
            coord27 = (self._upsCoord(coord26[1]), coord26[1])
            self.lithosSphericMantleShape = [
                coord23,
                coord24,
                coord25,
                coord26,
                coord27,
                self.coord5,
            ]
        elif self._isNodeOnEqualTriangle(coord26[1]) == False:
            print("3")
            coord27 = (self._upsCoord(coord26[1]), coord26[1])
            self.lithosSphericMantleShape = [
                coord23,
                coord24,
                coord25,
                coord26,
                coord27,
            ]

        # farBackarcLithosphericMantleShape
        self.lithoSphericMantleShapeFarBackarc = [coord24, coord22, coord28, coord25]
        print(coord23)

    def _isNodeOnEqualTriangle(self, modelHeigth: float) -> Tuple:
        if modelHeigth < self.coord5[1]:
            x = self._upsCoord(modelHeigth)
            return False
        else:
            x = self._equalLeggedTriangle(modelHeigth)
            return True

    def _upsCoord(self, y) -> float:
        y = y * self._magnitude
        x = (y - self.ups_b) / self.ups_a
        return x / self._magnitude

    def _equalLeggedTriangle(self, y) -> float:
        y = y * self._magnitude
        x = (y - self.elt_b) / self.elt_a
        return x / self._magnitude

    def getCrustSlabShapeArray(self) -> List[Tuple]:
        return self.crustSlabShape

    def getLithosphericMantleShapeForeArc(self) -> List[Tuple]:
        return self.lithosSphericMantleShape

    def getLithosphericMantleShapeFarBackArc(self) -> List[Tuple]:
        return self.lithoSphericMantleShapeFarBackarc

    def getUpperSlabShapeArray(self) -> List[Tuple]:
        return self.upperSlabPolygon

    def getMiddleSlabShapeArray(self) -> List[Tuple]:
        return self.middleSlabPolygon

    def getLowerSlabShapeArray(self) -> List[Tuple]:
        return self.lowerSlabPolygon
