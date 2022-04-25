import logging
import math
import os
import pickle
from time import time
from typing import Tuple

from underworld import conditions
from underworld import function as fn
from underworld import mesh, mpi, swarm, systems, utils
from underworld.function._function import Function
from underworld.scaling import units as u

from BaseModel import BaseModel
from CheckPointManager import CheckPointManager
from FigureManager import FigureManager
from modelParameters import ScalingCoefficientType
from modelParameters._Model_parameter_map import ModelParameterMap
from PlatePolygons import SubductionZonePolygons
from RheologyFunctions import RheologyFunctions


class SubductionModel(BaseModel):
    
    def __init__(self,
                 modelParameters: ModelParameterMap,
                 totalSteps: int,
                 checkPointSteps: int,
                 resolution:Tuple,
                 subductionZonePolygons:SubductionZonePolygons) -> None:
        super().__init__(modelParameters, totalSteps, checkPointSteps)
        self._resolution = resolution
        self._polygons = SubductionZonePolygons
        
    
    def _setMesh(self):
        self.mesh = mesh.FeMesh_Cartesian(
                elementType="Q1/dQ0",
                elementRes=(self.resolution),
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
    
    
    def _init_temperature_variables(self):
        self._temperatureField = mesh.MeshVariable(mesh=self.mesh, nodeDofCount=1)
        self._temperatureDotField = mesh.MeshVariable(mesh=self.mesh, nodeDofCount=1)
        self._temperatureDotField[...] = 0.
        
        proxyTemp = self.swarm.add_variable(dataType="double", count=1)
        proxyTemp.data[:] = 1.0

        for index in range(len(self.swarm.particleCoordinates.data)):
            coord = self.swarm.particleCoordinates.data[index][:]
            # if coord[1] < self.parameters.lowerMantleHeigth.nonDimensionalValue.magnitude:
            #     self.materialVariable.data[index] = self.lowerMantleIndex
            if self.slabUpperPoly.evaluate(tuple(coord)):
                proxyTemp.data[index] = 0.0
            if self.slabCorePoly.evaluate(tuple(coord)):
                proxyTemp.data[index] = 0.0
            if self.slabLowerPoly.evaluate(tuple(coord)):
                proxyTemp.data[index] = 0.0
        TmapSolver = utils.MeshVariable_Projection(self._temperatureField, proxyTemp)
        TmapSolver.solve()
        print("solved TemperatureField")
    
    @property
    def strainRate2ndInvariant(self) -> fn.Function:
        strainRateSecondInvariant = fn.tensor.second_invariant(
            fn.tensor.symmetric(self.velocityField.fn_gradient)
        )
        minimalStrainRate = fn.misc.constant(
            self.parameters.minimalStrainRate.nonDimensionalValue.magnitude
        )
        defaultStrainRate = fn.misc.constant(
            self.parameters.defaultStrainRate.nonDimensionalValue.magnitude
        )

        condition1 = [
            (self._solutionExists, strainRateSecondInvariant),
            (True, defaultStrainRate),
        ]

        existingStrainRate = fn.branching.conditional(condition1)

        condition2 = [
            (existingStrainRate <= minimalStrainRate, minimalStrainRate),
            (True, existingStrainRate),
        ]
        strainRate2ndInvariant = fn.branching.conditional(condition2)
        return strainRate2ndInvariant
    
    @property
    def vonMisesUpperLayerSP(self):
        sigmaY = (
            self.parameters.yieldStressOfSpTopLayer.nonDimensionalValue.magnitude
        )
        strainRateSecondInvariant = self.strainRate2ndInvariant
        effectiveViscosity = 0.5 * sigmaY / (strainRateSecondInvariant)
        return effectiveViscosity
    
    @property
    def stressFn(self) -> Function:
        return 2. * self.viscosityField * self.strainRateField
    
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
            [(self.depthFn > maxDepth, alteredViscosity), (self.depthFn < maxDepth, visTopLayer)]
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
        return fn.branching.map(
            fn_key=self.materialVariable, mapping=viscosityMap
        )
    
        
class OldSubductionModel:
    def __init__(
        self,
        *,
        name: str,
        modelParameterMap: ModelParameterMap,
        resolution: Tuple,
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
        self.resolution = resolution
        self.mesh = None
        # self.dissipation = self.swarm.add_variable(dataType="double", count=1)
        # self.storedEnergyRate = self.swarm.add_variable(dataType="double", count=1)

        self.totalSteps = totalSteps
        self.stepAmountCheckpoint = stepAmountCheckpoint
        self._setOutputPath()

        if fromCheckpoint:
            if fromCheckpointStep is None:
                raise ValueError

            self._initFromCheckPoint(fromCheckpointStep)
        else:
            if subductionZonePolygons is None:
                raise ValueError

            self.subductionZonePolygons = subductionZonePolygons
            self._initDefault()

    def _initDefault(self):

        self._setMesh()
        mpi.barrier()
        self._initSwarm()
        self.materialVariable = self.swarm.add_variable(dataType="int", count=1)
        self.previousStress = self.swarm.add_variable(dataType="double", count=3)
        self.viscosityField = self.swarm.add_variable(dataType="doubtle", count=1)

        self.previousStress.data[:] = [0.0, 0.0, 0.0]
        mpi.barrier()
        self._setupFields()

        self.rheologyCalculations = RheologyFunctions(self.parameters)
        mpi.barrier()
        self.meshHandle = None
        self._setupMaterialVarIndices()
        mpi.barrier()
        self._assignPolygons()
        mpi.barrier()
        self._fillTemperatureField()
        mpi.barrier()
        self._assignMaterialToVar()
        mpi.barrier()
        self._setBoundaryConditions()
        mpi.barrier()
        print("setBoundary")
        self._assignViscosityAndCreateMap()
        mpi.barrier()
        print("createdVis")
        self._assignStressAndCreateMap()
        mpi.barrier()
        print("createdStress")
        self._setBuoyancy()
        mpi.barrier()
        print("setBuoy")

        self._setAdvectionDiffusionSystem()
        mpi.barrier()
        print("advDi")
        self._setSwarmAdvectionSystem()
        mpi.barrier()
        print("swarmAd")
        self._setStokesSystem()
        mpi.barrier()
        print("setStokes")
        self._setStokesSolver()
        mpi.barrier

    def _getViscosityField(self):
        self.viscosityField.data[...] = self.viscosityFn.evaluate(self.swarm)

    def _initFromCheckPoint(self, step):
        self._setMesh()
        mpi.barrier()
        manager = CheckPointManager(self.name, self.outputPath)
        self.currentStep = step

        # self.mesh = manager.getMesh(self.preMesh)
        self.swarm = swarm.Swarm(mesh=self.mesh)
        manager.getSwarm(self.swarm, step)
        self.temperatureDotField = manager.getTemperatureDotField(step, self.mesh)
        self.materialVariable = manager.getMaterialVariable(step, self.swarm)
        self.previousStress = manager.getPreviousStress(step, self.swarm)
        self.pressureField = manager.getPressureField(step, self.mesh)
        self.velocityField = manager.getVelocityField(step, self.mesh)
        self.temperatureField = manager.getTemperatureField(step, self.mesh)
        self.meshHandle = None
        self.currentTime = manager.getLastTime(step)
        self.rheologyCalculations = RheologyFunctions(self.parameters)

        self._setupMaterialVarIndices()
        self._setBoundaryConditions()
        self._assignViscosityAndCreateMap()
        self._assignStressAndCreateMap()
        self._setBuoyancy()
        self._setAdvectionDiffusionSystem()
        self._setSwarmAdvectionSystem()
        self._setStokesSystem()
        self._setStokesSolver()

    def _setMesh(self):
        if self.mesh is None:
            

    def _initSwarm(self):

        

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
        self.lateralWalls = (
            self.mesh.specialSets["Top_VertexSet"]
            + self.mesh.specialSets["Bottom_VertexSet"]
        )
        self.rightWall = self.mesh.specialSets["Right_VertexSet"]
        self.leftWall = self.mesh.specialSets["Left_VertexSet"]
        self.topWall = self.mesh.specialSets["Top_VertexSet"]
        self.bottomWall = self.mesh.specialSets["Bottom_VertexSet"]

        self.VelocityBoundaryCondition = conditions.DirichletCondition(
            variable=self.velocityField,
            indexSetsPerDof=(self.verticalWalls, self.lateralWalls),
        )
        self.temperatureBoundaryCondition = conditions.DirichletCondition(
            variable=self.temperatureField,
            indexSetsPerDof=(self.topWall + self.verticalWalls),
        )

    def _setupFields(self):
        self.velocityField = mesh.MeshVariable(mesh=self.mesh, nodeDofCount=2)
        self.pressureField = mesh.MeshVariable(mesh=self.mesh.subMesh, nodeDofCount=1)
        self.temperatureField = mesh.MeshVariable(mesh=self.mesh, nodeDofCount=1)
        self.temperatureDotField = mesh.MeshVariable(mesh=self.mesh, nodeDofCount=1)
        self.strainRateField = mesh.MeshVariable(mesh=self.mesh.subMesh, nodeDofCount=1)
        self.strainRateField.data[...] = 0.0
        self.temperatureDotField.data[:] = 0.0
        self.velocityField.data[:] = 0.0
        self.pressureField.data[:] = 0.0

    def _setupMaterialVarIndices(self):
        self.upperMantleIndex = 0
        self.upperSlabIndex = 1
        self.lowerSlabIndex = 2
        # self.lowerMantleIndex = 4
        self.coreSlabIndex = 3

    def _assignPolygons(self):
        upperSlabShape = self.subductionZonePolygons.getUpperSlabShapeArray()
        lowerSlabShape = self.subductionZonePolygons.getLowerSlabShapeArray()
        middleSlabShape = self.subductionZonePolygons.getMiddleSlabShapeArray()
        self.slabUpperPoly = fn.shape.Polygon(upperSlabShape)
        self.slabLowerPoly = fn.shape.Polygon(lowerSlabShape)
        self.slabCorePoly = fn.shape.Polygon(middleSlabShape)

    def _fillTemperatureField(self):

        proxyTemp = self.swarm.add_variable(dataType="double", count=1)
        proxyTemp.data[:] = 1.0

        for index in range(len(self.swarm.particleCoordinates.data)):
            coord = self.swarm.particleCoordinates.data[index][:]
            # if coord[1] < self.parameters.lowerMantleHeigth.nonDimensionalValue.magnitude:
            #     self.materialVariable.data[index] = self.lowerMantleIndex
            if self.slabUpperPoly.evaluate(tuple(coord)):
                proxyTemp.data[index] = 0.0
            if self.slabCorePoly.evaluate(tuple(coord)):
                proxyTemp.data[index] = 0.0
            if self.slabLowerPoly.evaluate(tuple(coord)):
                proxyTemp.data[index] = 0.0
        TmapSolver = utils.MeshVariable_Projection(self.temperatureField, proxyTemp)
        TmapSolver.solve()
        print("solved TemperatureField")

    def _assignMaterialToVar(self):
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

    def _getDepthFunction(self):
        coordinate = fn.input()
        depthFn = self.mesh.maxCoord[1] - coordinate[1]
        return depthFn

    def _assignViscosityAndCreateMap(self):
        fnDepth = self._getDepthFunction()

        visTopLayer = fn.misc.min(
            self.rheologyCalculations.getEffectiveViscosityOfUpperLayerVonMises(
                self.velocityField
            ),
            self.parameters.spTopLayerViscosity.nonDimensionalValue.magnitude,  # self.parameters.spTopLayerViscosity.nonDimensionalValue.magnitude
        )

        maxDepth = 200e3 * u.meter
        maxDepth = self.parameters.scalingCoefficient.scalingForLength(
            maxDepth
        ).magnitude
        alteredViscosity = 50.0

        conditionTopLayer = fn.branching.conditional(
            [(fnDepth > maxDepth, alteredViscosity), (fnDepth < maxDepth, visTopLayer)]
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
        self.viscosityFn = fn.branching.map(
            fn_key=self.materialVariable, mapping=viscosityMap
        )

    def _assignStressAndCreateMap(self):
        Te = self.rheologyCalculations.getSymmetricStrainRateTensor(self.velocityField)
        visFn = self.viscosityFn

        viscousStressFn = 2.0 * visFn * Te

        stressMap = {
            self.upperMantleIndex: viscousStressFn,
            # self.lowerMantleIndex: viscousStressFn,
            self.lowerSlabIndex: viscousStressFn,
            self.coreSlabIndex: viscousStressFn,
            self.upperSlabIndex: viscousStressFn,
        }
        self.stressFn = fn.branching.map(
            fn_key=self.materialVariable, mapping=stressMap
        )
        self.stress2ndInvariant = fn.tensor.second_invariant(self.stressFn)

    def _setBuoyancy(self):
        print(self.temperatureField.data)
        ez = (0.0, -1.0)
        Ra = self.rheologyCalculations.getRayleighNumber()
        thermalDensityFn = Ra * (1.0 - self.temperatureField)
        self.buoyancyMapFn = thermalDensityFn * ez

    def _setStokesSystem(self):
        self.stokes = systems.Stokes(
            velocityField=self.velocityField,
            pressureField=self.pressureField,
            fn_bodyforce=self.buoyancyMapFn,
            fn_viscosity=self.viscosityFn,
            conditions=[
                self.VelocityBoundaryCondition,
            ],
        )

    def _setAdvectionDiffusionSystem(self):
        self.advectionDiffusion = systems.AdvectionDiffusion(
            phiField=self.temperatureField,
            phiDotField=self.temperatureDotField,
            velocityField=self.velocityField,
            fn_sourceTerm=0.0,
            fn_diffusivity=self.parameters.thermalDiffusivity.nonDimensionalValue.magnitude,
            conditions=[
                self.temperatureBoundaryCondition,
            ],
            # allow_non_q1=True,
            # method="SLCN",
        )

    def _setSwarmAdvectionSystem(self):
        self.swarmAdvector = systems.SwarmAdvector(
            self.velocityField, self.swarm, order=2
        )

    def _setStokesSolver(self):
        # try:
        self.solver = systems.Solver(self.stokes)
        self.solver.set_inner_method("mumps")

    # except RuntimeError:
    #     self.solver = systems.Solver(self.stokes)

    def _update(self, time, step):
        dt = self.parameters.deltaTime.nonDimensionalValue.magnitude

        # if dt > self.parameters.timeScaleStress.nonDimensionalValue.magnitude:
        #     dt = self.parameters.timeScaleStress.nonDimensionalValue.magnitude

        self.advectionDiffusion.integrate(dt)
        self.swarmAdvector.integrate(dt, update_owners=True)
        self.swarm.update_particle_owners()
        self.populationControl.repopulate()

        dt = dt * self.parameters.scalingCoefficient.timeCoefficient.magnitude
        self.rheologyCalculations.strainRateSolutionExists.value = True
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
            temperatureDotField=self.temperatureDotField,
            materialVariable=self.materialVariable,
            previousStress=self.previousStress,
            velocityField=self.velocityField,
            pressureField=self.pressureField,
            temperatureField=self.temperatureField,
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
        check_start_time = time()
        while self.currentStep < self.totalSteps:
            # self.solver.set_penalty(100)

            self.solver.solve(
                nonLinearIterate=True, nonLinearTolerance=0.1, print_stats=True
            )
            if self.currentStep == 1:
                self._assignViscosityAndCreateMap()
                mpi.barrier()
                self._assignStressAndCreateMap()
                mpi.barrier()
                # self._setBuoyancy()

                self._setAdvectionDiffusionSystem()
                mpi.barrier()
                self._setSwarmAdvectionSystem()
                mpi.barrier()
                self._setStokesSystem()
                mpi.barrier()
                self._setStokesSolver()

            if (
                self.currentStep % self.stepAmountCheckpoint == 0
                or self.currentStep == self.totalSteps - 1
            ):
                self._checkpoint(self.currentStep, self.currentTime)

                Vrms = math.sqrt(velSquared.evaluate()[0] / area.evaluate()[0])
                print(
                    f"{self.name = }, {self.currentStep = } {self.currentTime = :.3e} {Vrms = :.3e} "
                )

            mpi.barrier()
            newTime, newStep = self._update(self.currentTime, self.currentStep)
            self.currentStep = newStep
            self.currentTime = newTime
            self.figureManager.incrementStoreStep()
            check_endTime = time()
            time_for_loop = check_endTime - check_start_time
            check_start_time = check_endTime
            print(f"{self.currentStep = }, {time_for_loop = }")

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
