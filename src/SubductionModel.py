import logging
import math
import os

from underworld import conditions
from underworld import function as fn
from underworld import mesh, mpi, swarm, systems, utils
from underworld.function._function import Function

from CheckPointManager import CheckPointManager
from FigureManager import FigureManager
from modelParameters import ScalingCoefficientType
from modelParameters._Model_parameter_map import ModelParameterMap
from PlatePolygons import SubductionZonePolygons
from RheologyFunctions import RheologyFunctions


class SubductionModel:
    def __init__(
        self,
        *,
        name: str,
        modelParameterMap: ModelParameterMap,
        totalSteps,
        stepAmountCheckpoint,
        subductionZonePolygons: SubductionZonePolygons = None,
        fromCheckpoint: bool = False,
        fromCheckpointStep=None,
    ) -> None:
        """
        If you want to continue from a checkpoint param: subducionZonePolygons can be None
        """
        self.name = name
        self.parameters = modelParameterMap
        self.currentStep = 0
        self.currentTime = 0.0

        # self.dissipation = self.swarm.add_variable(dataType="double", count=1)
        # self.storedEnergyRate = self.swarm.add_variable(dataType="double", count=1)

        self.totalSteps = totalSteps
        self.stepAmountCheckpoint = stepAmountCheckpoint
        self._setOutputPath()
        if fromCheckpoint:
            if fromCheckpointStep is None:
                raise ValueError

            self.currentStep = fromCheckpointStep
            self._initFromCheckPoint(self.currentStep)
        else:
            if subductionZonePolygons is None:
                raise ValueError

            self.subductionZonePolygons = subductionZonePolygons
            self._initDefault()

    def _initDefault(self):
        self._initMeshAndField()
        self.materialVariable = self.swarm.add_variable(dataType="int", count=1)
        self.previousStress = self.swarm.add_variable(dataType="double", count=3)
        self.previousStress.data[:] = [0.0, 0.0, 0.0]
        self._initializeFields()
        self.rheologyCalculations = RheologyFunctions(self.parameters)
        self.meshHandle = None
        self._initializeMaterialVariableDataIndices()
        self._initializePolygons()
        self._assignMaterialToParticles()
        self._setBoundaryConditions()
        self._initializeViscosityMapFn()
        self._initializeStressMapFn()
        self._setBuoyancy()
        self._initializeStokes()
        self._initializeAdvectionDiffusion()
        self._initializeSwarmAdvector()
        self._initializeVeStressMap()
        self._initializeSolver()

    def _initFromCheckPoint(self, step):
        manager = CheckPointManager(self.name, self.outputPath)

        self.mesh = manager.getMesh()
        self.swarm = manager.getSwarm(step)

        self.materialVariable = manager.getMaterialVariable(step, self.swarm)
        self.previousStress = manager.getPreviousStress(step, self.swarm)
        self.pressureField = manager.getPressureField(step, self.mesh)
        self.velocityField = manager.getVelocityField(step, self.mesh)
        self.temperatureField = manager.getTemperatureField(step, self.mesh)
        self.temperatureDotField = manager.getTemperatureDotField(step, self.mesh)
        self.meshHandle = None
        self.currentTime = manager.getLastTime(step)
        self.rheologyCalculations = RheologyFunctions(self.parameters)

        self._initializeMaterialVariableDataIndices()
        self._setBoundaryConditions()
        self._initializeViscosityMapFn()
        self._initializeStressMapFn()
        self._setBuoyancy()
        self._initializeStokes()
        self._initializeAdvectionDiffusion()
        self._initializeSwarmAdvector()
        self._initializeVeStressMap()
        self._initializeSolver()
        if mpi.rank == 0:
            self.figureManager.store.step = self.currentStep

    def _initMeshAndField(self):
        self.mesh = mesh.FeMesh_Cartesian(
            elementType="Q1/dQ0",
            elementRes=(200, 100),
            minCoord=(0.0, 0.0),
            maxCoord=(
                self.parameters.modelLength.nonDimensionalValue.magnitude,
                self.parameters.modelHeight.nonDimensionalValue.magnitude,
            ),
            periodic=[True, False],
        )
        self.swarm = swarm.Swarm(mesh=self.mesh)
        self.swarmLayout = swarm.layouts.PerCellSpaceFillerLayout(
            swarm=self.swarm, particlesPerCell=20
        )
        self.swarm.populate_using_layout(self.swarmLayout)

    def _setOutputPath(self):
        if mpi.rank == 0:
            try:
                os.mkdir("./output")

            except FileExistsError:
                pass
            try:
                os.mkdir(f"./output/{self.name}")
            except FileExistsError:
                pass
        self.outputPath = f"./output/{self.name}"
        self.figureManager = FigureManager(self.outputPath, self.name)

    def _setBoundaryConditions(self):
        self.verticalWalls = (
            self.mesh.specialSets["Left_VertexSet"]
            + self.mesh.specialSets["Right_VertexSet"]
        )
        self.horizontalWalls = (
            self.mesh.specialSets["Top_VertexSet"]
            + self.mesh.specialSets["Bottom_VertexSet"]
        )
        self.rightWall = self.mesh.specialSets["Right_VertexSet"]
        self.leftWall = self.mesh.specialSets["Left_VertexSet"]
        self.topWall = self.mesh.specialSets["Top_VertexSet"]
        self.bottomWall = self.mesh.specialSets["Bottom_VertexSet"]

        fixed = self.mesh.specialSets["Empty"]
        fixNode = self.rightWall + self.bottomWall
        fixed += fixNode

        self.VelocityBoundaryCondition = conditions.DirichletCondition(
            variable=self.velocityField,
            indexSetsPerDof=(self.verticalWalls, self.horizontalWalls),
        )
        self.temperatureBoundaryCondition = conditions.DirichletCondition(
            variable=self.temperatureField,
            indexSetsPerDof=(self.topWall + self.verticalWalls),
        )

    def _initializeFields(self):
        self.velocityField = mesh.MeshVariable(mesh=self.mesh, nodeDofCount=2)
        self.pressureField = mesh.MeshVariable(mesh=self.mesh.subMesh, nodeDofCount=1)
        self.temperatureField = mesh.MeshVariable(mesh=self.mesh, nodeDofCount=1)
        self.temperatureDotField = mesh.MeshVariable(mesh=self.mesh, nodeDofCount=1)
        self.velocityField.data[:] = [0.0, 0.0]
        self.pressureField.data[:] = 0.0
        self.temperatureDotField.data[:] = 0.0
        self.temperatureField.data[:] = 0

    def _initializeMaterialVariableDataIndices(self):
        self.upperMantleIndex = 0
        self.upperSlabIndex = 1
        self.lowerSlabIndex = 2
        # self.lowerMantleIndex = 4
        self.coreSlabIndex = 3

    def _initializePolygons(self):
        upperSlabShape = self.subductionZonePolygons.getUpperSlabShapeArray()
        lowerSlabShape = self.subductionZonePolygons.getLowerSlabShapeArray()
        middleSlabShape = self.subductionZonePolygons.getMiddleSlabShapeArray()

        self.slabUpperPoly = fn.shape.Polygon(upperSlabShape)
        self.slabLowerPoly = fn.shape.Polygon(lowerSlabShape)
        self.slabCorePoly = fn.shape.Polygon(middleSlabShape)

    def _assignMaterialToParticles(self):
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

    def _initializeViscosityMapFn(self):
        visTopLayer = (
            self.rheologyCalculations.getEffectiveViscosityOfUpperLayerVonMises(
                self.velocityField
            )
            # self.parameters.spTopLayerViscosity.nonDimensionalValue.magnitude
        )
        visCoreLayer = (
            self.rheologyCalculations.getEffectiveViscosityOfViscoElasticCore()
            # self.parameters.spCoreLayerViscosity.nonDimensionalValue.magnitude
        )
        visBottomLayer = (
            self.parameters.spBottomLayerViscosity.nonDimensionalValue.magnitude
        )
        toplayerViscosity = (
            self.parameters.spTopLayerViscosity.nonDimensionalValue.magnitude
        )

        afterStep = fn.exception.SafeMaths(
            fn.misc.min(
                visTopLayer,
                toplayerViscosity,
            )
        )
        conditionalTopLayer = fn.branching.conditional(
            (
                (self.currentStep >= 1, afterStep),
                (self.currentStep < 1, toplayerViscosity),
            )
        )

        # print(f"{var = }")
        viscosityMap = {
            self.upperMantleIndex: self.parameters.upperMantleViscosity.nonDimensionalValue.magnitude,
            # self.lowerMantleIndex: self.parameters.lowerMantleViscosity.nonDimensionalValue.magnitude,
            self.lowerSlabIndex: round(visBottomLayer, 1),
            self.coreSlabIndex: visCoreLayer,
            self.upperSlabIndex: conditionalTopLayer,
        }
        self.viscosityFn = fn.branching.map(
            fn_key=self.materialVariable, mapping=viscosityMap
        )

    def _initializeStressMapFn(self):
        Te = self.rheologyCalculations.getSymmetricStrainRateTensor(self.velocityField)
        visFn = self.viscosityFn

        viscousStressFn = 2.0 * visFn * Te

        visEff = self.rheologyCalculations.getEffectiveViscosityOfViscoElasticCore()
        shearMod = self.parameters.coreShearModulus.nonDimensionalValue.magnitude
        dt = self.parameters.timeScaleStress.nonDimensionalValue.magnitude

        self.elasticStressFn = visEff / (shearMod * dt) * self.previousStress

        self.viscoElasticStressFn = viscousStressFn + self.elasticStressFn

        stressMap = {
            self.upperMantleIndex: viscousStressFn,
            # self.lowerMantleIndex: viscousStressFn,
            self.lowerSlabIndex: viscousStressFn,
            self.coreSlabIndex: self.viscoElasticStressFn,
            self.upperSlabIndex: viscousStressFn,
        }
        self.stressFn = fn.branching.map(
            fn_key=self.materialVariable, mapping=stressMap
        )
        self.stress2ndInvariant = fn.tensor.second_invariant(self.stressFn)

    def _setBuoyancy(self):
        ez = (0.0, -1.0)
        Ra = self.rheologyCalculations.getRayleighNumber()
        thermalDensityFn = Ra * (1.0 - self.temperatureField)
        self.buoyancyMapFn = thermalDensityFn * ez

    def _initializeStokes(self):
        self.stokes = systems.Stokes(
            velocityField=self.velocityField,
            pressureField=self.pressureField,
            fn_bodyforce=self.buoyancyMapFn,
            fn_viscosity=self.viscosityFn,
            fn_stresshistory=self.elasticStressFn,
            conditions=[
                self.VelocityBoundaryCondition,
            ],
        )

    def _initializeAdvectionDiffusion(self):
        self.advectionDiffusion = systems.AdvectionDiffusion(
            phiField=self.temperatureField,
            phiDotField=self.temperatureDotField,
            velocityField=self.velocityField,
            fn_sourceTerm=0.0,
            fn_diffusivity=1.0,
            conditions=[
                self.temperatureBoundaryCondition,
            ],
            allow_non_q1=True,
        )

    def _initializeSwarmAdvector(self):
        self.swarmAdvector = systems.SwarmAdvector(
            self.velocityField, self.swarm, order=2
        )

    def _initializeVeStressMap(self):
        self.veStressMap = {
            self.upperMantleIndex: [0.0, 0.0, 0.0],
            self.upperSlabIndex: [0.0, 0.0, 0.0],
            self.lowerSlabIndex: [0.0, 0.0, 0.0],
            self.coreSlabIndex: self.viscoElasticStressFn,
        }
        self.veStressFn = fn.branching.map(
            fn_key=self.materialVariable, mapping=self.veStressMap
        )

    def _initializeSolver(self):
        # try:
        self.solver = systems.Solver(self.stokes)
        self.solver.set_inner_method("mumps")

    # except RuntimeError:
    #     self.solver = systems.Solver(self.stokes)

    def _update(self, time, step):
        dt = self.advectionDiffusion.get_max_dt()

        dt_e = self.parameters.timeScaleStress.nonDimensionalValue.magnitude

        if dt > (dt_e / 3.0):
            dt = dt_e / 3

        phi = dt / dt_e

        veStressFn_data = self.veStressFn.evaluate(self.swarm)

        self.previousStress.data[:] = (
            phi * veStressFn_data[:] + (1 - phi) * self.previousStress.data[:]
        )

        self.advectionDiffusion.integrate(dt)
        self.swarmAdvector.integrate(dt)
        dt = dt * self.parameters.scalingCoefficient.timeCoefficient.magnitude
        return time + dt, step + 1

    def getMeshHandle(self):
        if self.meshHandle is None:
            try:
                self.meshHandle = self.mesh.save(self.outputPath + "mesh.00000.h5")
            except FileExistsError:
                if mpi.rank == 0:
                    os.remove(self.outputPath + "mesh.00000.h5")
                mpi.barrier()
                return self.getMeshHandle()

        return self.meshHandle

    def _checkpoint(self, step, time):
        manager = CheckPointManager(self.name, self.outputPath)
        manager.checkPoint(
            step=step,
            swarm=self.swarm,
            mesh=self.mesh,
            materialVariable=self.materialVariable,
            previousStress=self.previousStress,
            velocityField=self.velocityField,
            pressureField=self.pressureField,
            temperatureField=self.temperatureField,
            temperatureDotField=self.temperatureDotField,
            figureManager=self.figureManager,
            meshHandle=self.getMeshHandle(),
            strainRate2ndInvariant=self.rheologyCalculations.getStrainRateSecondInvariant(
                self.velocityField
            ),
            viscosityFn=self.viscosityFn,
            stress2ndInvariant=self.stress2ndInvariant,
            time=time,
        )

    def run(self):
        # try:

        velSquared = utils.Integral(
            fn.math.dot(self.velocityField, self.velocityField), self.mesh
        )
        area = utils.Integral(1.0, self.mesh)

        while self.currentStep < self.totalSteps:
            # self.solver.set_penalty(100)
            self.solver.solve(nonLinearIterate=True, nonLinearTolerance=0.1)

            if (
                self.currentStep % self.stepAmountCheckpoint == 0
                or self.currentStep == self.totalSteps - 1
            ):
                self._checkpoint(self.currentStep, self.currentTime)

                Vrms = math.sqrt(velSquared.evaluate()[0] / area.evaluate()[0])
                logging.debug(
                    f"{self.name = }, {self.currentStep = } {self.currentTime = :.3e} {Vrms = :.3e} "
                )
            newTime, newStep = self._update(self.currentTime, self.currentStep)
            self.currentStep = newStep
            self.currentTime = newTime
            self.figureManager.incrementStoreStep()

    # except KeyboardInterrupt:
    # try:
    #     Vrms = math.sqrt(velSquared.evaluate()[0] / area.evaluate()[0])
    #     self._checkpoint(self.currentStep, self.currentTime)
    #     logging.debug(
    #         f"sucessfully checkpointed after keyboardInterrupt {self.currentStep =}, {self.currentTime = }, {Vrms = }"
    #     )
    # except Exception as e:
    #     logging.exception(
    #         f" Failed final checkpoint after KeyboardInterrupt {e = }"
    #     )
    #     raise e

    # except Exception as e:
    #     logging.exception(f" Run Failed {e = }", stack_info=True)
    # try:
    #     Vrms = math.sqrt(velSquared.evaluate()[0] / area.evaluate()[0])
    #     self._checkpoint(self.currentStep, self.currentTime)
    #     logging.debug(
    #         f"sucessfully checkpointed after exception {self.currentStep =}, {self.currentTime = }, {Vrms = }"
    #     )
    #     raise e
    # except Exception as a:
    #     logging.exception(f" Failed final checkpoint {a = }")
    #     raise a


# TODO create a start from checkpoint function
