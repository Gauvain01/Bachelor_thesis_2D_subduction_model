import math
from time import time
from typing import Tuple
from unicodedata import name

from underworld import conditions
from underworld import function as fn
from underworld import mesh, mpi, swarm, utils
from underworld.function._function import Function
from underworld.scaling import units as u

from BaseModel import BaseModel
from modelParameters._Model_parameter_map import ModelParameterMap
from PlatePolygons import SubductionZonePolygons


class SubductionModel(BaseModel):
    def __init__(
        self,
        modelParameters: ModelParameterMap,
        totalSteps: int,
        checkPointSteps: int,
        resolution: Tuple,
        subductionZonePolygons: SubductionZonePolygons,
        name: str,
        restart: bool = False,
        restartStep: int = None,
    ) -> None:
        self._resolution = resolution
        self._polygons = subductionZonePolygons
        super().__init__(
            modelParameters, totalSteps, checkPointSteps, name, restart, restartStep
        )

    def _setMesh(self):
        self.mesh = mesh.FeMesh_Cartesian(
            elementType="Q1/dQ0",
            elementRes=(self._resolution),
            minCoord=(0.0, 0.0),
            maxCoord=(
                self.parameters.modelLength.nonDimensionalValue.magnitude,
                self.parameters.modelHeight.nonDimensionalValue.magnitude,
            ),
        )

    def _setSwarm(self):
        self.swarm = swarm.Swarm(mesh=self.mesh, particleEscape=True)
        self.swarm.allow_parallel_nn = True
        self.swarmLayout = swarm.layouts.PerCellSpaceFillerLayout(
            swarm=self.swarm, particlesPerCell=20
        )
        self.swarm.populate_using_layout(self.swarmLayout)
        self.populationControl = swarm.PopulationControl(
            self.swarm,
            particlesPerCell=20,
            aggressive=True,
            splitThreshold=0.15,
            maxSplits=10,
        )

    def _initTemperatureVariables(self):

        self._temperatureDotField.data[...] = 0.0

        self._proxyTemp.data[...] = 1.0

        for index in range(len(self.swarm.particleCoordinates.data)):
            coord = self.swarm.particleCoordinates.data[index][:]
            # if coord[1] < self.parameters.lowerMantleHeigth.nonDimensionalValue.magnitude:
            #     self.materialVariable.data[index] = self.lowerMantleIndex
            if self.slabUpperPoly.evaluate(tuple(coord)):
                self._proxyTemp.data[index] = 0.0
            if self.slabCorePoly.evaluate(tuple(coord)):
                self._proxyTemp.data[index] = 0.0
            if self.slabLowerPoly.evaluate(tuple(coord)):
                self._proxyTemp.data[index] = 0.0
        TmapSolver = utils.MeshVariable_Projection(
            self._temperatureField, self._proxyTemp
        )
        TmapSolver.solve()
        print("solved TemperatureField")

    @property
    def slabUpperPoly(self):
        upper = self._polygons.getUpperSlabShapeArray()
        return fn.shape.Polygon(upper)

    @property
    def slabCorePoly(self):
        core = self._polygons.getMiddleSlabShapeArray()
        return fn.shape.Polygon(core)

    @property
    def slabLowerPoly(self):
        lower = self._polygons.getLowerSlabShapeArray()
        return fn.shape.Polygon(lower)

    @property
    def strainRate2ndInvariant(self) -> fn.Function:
        strainRateSecondInvariant = fn.tensor.second_invariant(
            fn.tensor.symmetric(self.velocityField.fn_gradient)
        )
        minimalStrainRate = (
            self.parameters.minimalStrainRate.nonDimensionalValue.magnitude
        )

        defaultStrainRate = (
            self.parameters.defaultStrainRate.nonDimensionalValue.magnitude
        )

        conditionToCheckForExistingSolve = [
            (self._solutionExists, strainRateSecondInvariant),
            (True, defaultStrainRate),
        ]

        existingStrainRate = fn.branching.conditional(conditionToCheckForExistingSolve)

        conditionForMinimalStrain = [
            (existingStrainRate <= minimalStrainRate, minimalStrainRate),
            (True, existingStrainRate),
        ]
        strainRate2ndInvariant = fn.branching.conditional(conditionForMinimalStrain)
        return strainRate2ndInvariant

    @property
    def vonMisesUpperLayerSP(self):
        sigmaY = self.parameters.yieldStressOfSpTopLayer.nonDimensionalValue.magnitude
        strainRateSecondInvariant = self.strainRate2ndInvariant
        effectiveViscosity = 0.5 * sigmaY / strainRateSecondInvariant
        return effectiveViscosity

    @property
    def stressFn(self) -> Function:
        return 2.0 * self.viscosityField * self.strainRateField

    @property
    def depthFn(self) -> Function:
        coordinate = fn.input()
        depthFn = self.mesh.maxCoord[1] - coordinate[1]
        return depthFn

    @property
    def viscosityFn(self) -> Function:

        visTopLayer = fn.misc.min(
            self.vonMisesUpperLayerSP,
            self.parameters.spTopLayerViscosity.nonDimensionalValue.magnitude,  # self.parameters.spTopLayerViscosity.nonDimensionalValue.magnitude
        )

        maxDepth = 200e3 * u.meter
        maxDepth = self.parameters.scalingCoefficient.scalingForLength(
            maxDepth
        ).magnitude
        alteredViscosity = 50.0

        conditionTopLayer = fn.branching.conditional(
            [
                (self.depthFn > maxDepth, alteredViscosity),
                (self.depthFn < maxDepth, visTopLayer),
            ]
        )

        visCoreLayer = (
            self.parameters.spCoreLayerViscosity.nonDimensionalValue.magnitude
            # self.parameters.spCoreLayerViscosity.nonDimensionalValue.magnitude
        )
        visBottomLayer = (
            self.parameters.spBottomLayerViscosity.nonDimensionalValue.magnitude
        )

        # print(f"{var = }")
        viscosityMap = {
            self.upperMantleIndex: self.parameters.upperMantleViscosity.nonDimensionalValue.magnitude,
            # self.lowerMantleIndex: self.parameters.lowerMantleViscosity.nonDimensionalValue.magnitude,
            self.lowerSlabIndex: round(visBottomLayer, 1),
            self.coreSlabIndex: visCoreLayer,
            self.upperSlabIndex: conditionTopLayer,
        }
        return fn.branching.map(fn_key=self.materialVariable, mapping=viscosityMap)

    @property
    def rayleighNumber(self):

        ls = self.parameters.modelHeight.dimensionalValue.magnitude
        rhoRef = self.parameters.referenceDensity.dimensionalValue.magnitude

        g = self.parameters.gravitationalAcceleration.dimensionalValue.magnitude
        alpha = self.parameters.thermalExpansivity.dimensionalValue.magnitude

        deltaT = self.parameters.temperatureContrast.dimensionalValue.magnitude

        k = self.parameters.thermalDiffusivity.dimensionalValue.magnitude

        eta = self.parameters.referenceViscosity.dimensionalValue.magnitude

        rayLeighNumber = ((ls**3) * rhoRef * g * alpha * deltaT) / (k * eta)
        return rayLeighNumber

    @property
    def buoyancyFn(self) -> Function:
        ez = (0.0, -1.0)
        Ra = self.rayleighNumber
        thermalDensityFn = Ra * (1.0 - self.temperature)
        buoyancyMapFn = thermalDensityFn * ez
        return buoyancyMapFn

    @property
    def velocityBC(self):
        verticalWalls = (
            self.mesh.specialSets["Left_VertexSet"]
            + self.mesh.specialSets["Right_VertexSet"]
        )
        lateralWalls = (
            self.mesh.specialSets["Top_VertexSet"]
            + self.mesh.specialSets["Bottom_VertexSet"]
        )

        VelocityBoundaryCondition = conditions.DirichletCondition(
            variable=self.velocityField,
            indexSetsPerDof=(verticalWalls, lateralWalls),
        )
        return VelocityBoundaryCondition

    @property
    def temperatureBC(self):
        verticalWalls = (
            self.mesh.specialSets["Left_VertexSet"]
            + self.mesh.specialSets["Right_VertexSet"]
        )
        topWall = self.mesh.specialSets["Top_VertexSet"]
        temperatureBoundaryCondition = conditions.DirichletCondition(
            variable=self.temperature,
            indexSetsPerDof=(topWall + verticalWalls),
        )
        return temperatureBoundaryCondition

    @property
    def outputPath(self):
        path = f"./output/{self.name}"
        return path

    @property
    def upperMantleIndex(self):
        return 0

    @property
    def upperSlabIndex(self):
        return 1

    @property
    def coreSlabIndex(self):
        return 2

    @property
    def lowerSlabIndex(self):
        return 3

    def _initMaterialVariable(self):
        self.materialVariable.data[:] = self.upperMantleIndex

        for index in range(len(self.swarm.particleCoordinates.data)):
            coord = self.swarm.particleCoordinates.data[index][:]
            # if coord[1] < self.parameters.lowerMantleHeigth.nonDimensionalValue.magnitude:
            #     self.materialVariable.data[index] = self.lowerMantleIndex
            if self.slabUpperPoly.evaluate(tuple(coord)):
                self.materialVariable.data[index] = self.upperSlabIndex
            if self.slabCorePoly.evaluate(tuple(coord)):
                self.materialVariable.data[index] = self.coreSlabIndex
            if self.slabLowerPoly.evaluate(tuple(coord)):
                self.materialVariable.data[index] = self.lowerSlabIndex

    def _update(self):
        dt = self.parameters.deltaTime.nonDimensionalValue.magnitude

        self.advectionDiffusionSystem.integrate(dt)
        self.swarmAdvector.integrate(dt, update_owners=True)
        self.swarm.update_particle_owners()
        self.populationControl.repopulate()

        self._solutionExists.value = True
        self.modelTime += dt
        self.modelStep += 1
        self.figureManager.incrementStoreStep()

    @property
    def vrms(self):
        velSquared = utils.Integral(
            fn.math.dot(self.velocityField, self.velocityField), self.mesh
        )
        area = utils.Integral(1.0, self.mesh)
        vrms = math.sqrt(velSquared.evaluate()[0] / area.evaluate()[0])
        return vrms

    def run(self):

        check_start_time = time()
        while self.modelStep < self.totalSteps:

            self.solver.solve(
                nonLinearIterate=True, nonLinearTolerance=0.1, print_stats=True
            )

            if (
                self.modelStep % self.checkPointSteps == 0
                or self.modelStep == self.totalSteps - 1
            ):
                self._checkPoint()

                print(
                    f"{self.name = }, {self.modelStep = } {self.modelTime = :.3e} {self.vrms = :.3e} "
                )

            mpi.barrier()
            self._update()
            check_endTime = time()
            time_for_loop = check_endTime - check_start_time
            check_start_time = check_endTime
            print(f"{self.modelStep = }, {time_for_loop = }")
