from __future__ import annotations

import json
import os
from abc import abstractmethod

from underworld import function as fn
from underworld import mesh, mpi
from underworld import scaling as sca
from underworld import swarm, systems
from underworld.function import Function

from CheckPointManager import CheckPointManager
from FigureManager import FigureManager
from meshProjector import meshProjector
from modelParameters import ModelParameterMap


class BaseModel:
    def __init__(
        self,
        modelParameters: ModelParameterMap,
        totalSteps: int,
        checkPointSteps: int,
        name: str,
        restart: bool = False,
        restartStep: int = None,
    ) -> None:
        self.restart = restart
        self.restartStep = restartStep
        self.name = name
        self.totalSteps = totalSteps
        self.checkPointSteps = checkPointSteps
        self.parameters = modelParameters
        self.modelStep = 0
        self.modelTime = 0
        self._swarmVarForSaving = []
        self._meshVarForSaving = []
        self._solutionExists = fn.misc.constant(False)
        self._hasTemperatureDotField = False
        self._hasTemperatureField = False
        self._figureManager = FigureManager(self.outputPath, self.name)
        self._setMesh()
        self._setSwarm()
        mpi.barrier()
        self.addMeshVariable(
            "pressureField", "double", 1, subMesh=True, restartVariable=True
        )
        self.addSwarmVariable(
            "_materialVariable", dataType="int", count=1, restartVariable=True
        )
        self.addMeshVariable("velocityField", "double", 2, restartVariable=True)
        self.addMeshVariable("_strainRateField", "double", 1, subMesh=True)
        self.addSwarmVariable("_viscosityField", "double", 1)
        self.addMeshVariable("_projectedViscosity", "double", 1, subMesh=True)
        self.addSwarmVariable("_stressField", "double", 1, restartVariable=True)

        self.pressureField.data[:] = 0.0
        self.addSwarmVariable("_proxyTemp", "double", 1)
        self.addMeshVariable(
            "_temperatureField", "double", nodeDofCount=1, restartVariable=True
        )
        self.addMeshVariable("_temperatureDotField", "double", nodeDofCount=1)
        self._initMaterialVariable()
        self._initTemperatureVariables()
        mpi.barrier()
        self.testStokes()
        self.addMeshVariable(
            "_projectedStressField", "double", nodeDofCount=1, subMesh=True
        )
        self.addSwarmVariable("_stressTensor", "double", count=3)
        self.addMeshVariable(
            "_projectedStressTensor", "double", nodeDofCount=3, subMesh=True
        )

        self._makeOutputDir()
        self.testStokes()

    def _makeOutputDir(self):
        if mpi.rank == 0:
            try:
                os.mkdir(self.outputPath)
            except FileExistsError:
                pass
        mpi.barrier()

    def addMeshVariable(
        self,
        name: str,
        dataType: str,
        nodeDofCount,
        subMesh: bool = False,
        restartVariable: bool = False,
    ):
        if subMesh:
            field = mesh.MeshVariable(self.mesh.subMesh, nodeDofCount, dataType)
        else:
            field = mesh.MeshVariable(self.mesh, nodeDofCount, dataType)
        if restartVariable:
            self._meshVarForSaving.append((name, subMesh))
        if self.restart and restartVariable:
            self.checkPointManager.loadField(
                name=name, step=self.restartStep, field=field
            )
        setattr(self, name, field)
        self._meshVarForSaving.append((name, subMesh))

    def addSwarmVariable(
        self, name: str, dataType: str, count: int, restartVariable: bool = False
    ):
        self.swarm: swarm.Swarm
        field = self.swarm.add_variable(dataType, count)
        if restartVariable:
            self._swarmVarForSaving.append(name)
        if self.restart and restartVariable:
            self.checkPointManager.loadField(
                name=name, step=self.restartStep, field=field
            )
        setattr(self, name, field)

    @abstractmethod
    def _initTemperatureVariables(self):
        pass

    @abstractmethod
    def _setMesh(self):
        pass

    @abstractmethod
    def _setSwarm(self):
        pass

    @property
    @abstractmethod
    def outputPath(self):
        pass

    @abstractmethod
    def _initMaterialVariable(self):
        pass

    @property
    @abstractmethod
    def stressFn(self) -> Function:
        pass

    @property
    @abstractmethod
    def strainRate2ndInvariant(self):
        pass

    @property
    @abstractmethod
    def viscosityFn(self) -> Function:
        pass

    @property
    @abstractmethod
    def buoyancyFn(self) -> Function:
        pass

    @property
    @abstractmethod
    def velocityBC(self):
        pass

    @property
    @abstractmethod
    def temperatureBC(self):
        pass

    @property
    def materialVariable(self):
        return self._materialVariable

    @property
    def checkPointManager(self) -> CheckPointManager:
        obj = CheckPointManager(self.name, self.outputPath)
        return obj

    @property
    def temperature(self):
        if not self._hasTemperatureField:
            if hasattr(self, "_temperatureField"):
                self._hasTemperatureField = True
            else:
                raise AttributeError("_temperatureField not assigned")
        return self._temperatureField

    @property
    def temperatureDotField(self):
        if not self._hasTemperatureDotField:
            if hasattr(self, "_temperatureDotField"):
                self._hasTemperatureDotField = True
            else:
                raise AttributeError("_temperatureDotField not assigned")

        return self._temperatureDotField

    @property
    def viscosityField(self):
        self._viscosityField.data[...] = self.viscosityFn.evaluate(self.swarm)
        return self._viscosityField

    @property
    def projectedViscosityField(self):
        meshProjector(self._projectedViscosity, self.viscosityField)
        return self._projectedViscosity

    @property
    def strainRateField(self):
        self._strainRateField.data[:] = self.strainRate2ndInvariant.evaluate(
            self.mesh.subMesh
        )
        return self._strainRateField

    @property
    def projectedStressField(self):

        self._stressField.data[...] = self.stressFn.evaluate(self.swarm)
        meshProjector(self._projectedStressField, self._stressField)
        return self._projectedStressField

    @property
    def projectedStressTensor(self):
        self._stressTensor.data[...] = self.stressFn.evaluate(self.swarm)
        meshProjector(self._projectedStressTensor, self._stressTensor)
        return self._projectedStressTensor

    def testStokes(self):

        velField = self.velocityField
        pField = self.pressureField

        viscosityFn = self.viscosityFn

        bodyForce = self.buoyancyFn

        condition = self.velocityBC

        print(" i am here")
        mpi.barrier()
        systems.Stokes(
            velocityField=velField,
            pressureField=pField,
            fn_viscosity=viscosityFn,
            fn_bodyforce=bodyForce,
            conditions=[
                condition,
            ],
        )
        print("tested stokes succesfully")

    @property
    def stokes(self):

        self._stokes = systems.Stokes(
            velocityField=self.velocityField,
            pressureField=self.pressureField,
            fn_viscosity=self.viscosityFn,
            fn_bodyforce=self.buoyancyFn,
            conditions=[
                self.velocityBC,
            ],
        )
        print("i got the stokes")
        return self._stokes

    @property
    def swarmAdvector(self):
        obj = systems.SwarmAdvector(self.velocityField, self.swarm, order=2)
        return obj

    @property
    def solver(self):
        mpi.barrier()
        print(" i past the barrier for the solver")
        stokes = self.stokes
        self._solver = systems.Solver(stokes)
        self._solver.set_inner_method("mumps")
        print("i got the solver")
        return self._solver

    @property
    def meshHandle(self):
        with mpi.call_pattern("collective"):
            try:
                handle = self._meshHandle
            except AttributeError:
                self._meshHandle = self.mesh.save(self.outputPath + "mesh.00000.h5")
                handle = self._meshHandle
            except Exception as e:
                raise ValueError(f"problem with meshHandle {e}")

            return handle

    @property
    def swarmHandle(self):
        with mpi.call_pattern("collective"):
            try:
                handle = self._swarmHandle
            except AttributeError:
                self._swarmHandle = self.swarm.save(self.outputPath + "mesh.00000.h5")
                handle = self._swarmHandle
            except Exception as e:
                raise ValueError(f"problem with swarmHandle {e}")

            return handle

    @property
    def advectionDiffusionSystem(self):
        obj = systems.AdvectionDiffusion(
            phiField=self.temperature,
            fn_diffusivity=0.0,
            fn_sourceTerm=0.0,
            phiDotField=self.temperatureDotField,
            conditions=[
                self.temperatureBC,
            ],
            velocityField=self.velocityField,
        )

        return obj

    @property
    def figureManager(self):
        return self._figureManager

    def solve(self):
        self.solver.solve(
            nonLinearIterate=True, nonLinearTolerance=1e-2, print_stats=True
        )

    @abstractmethod
    def _update(self):
        pass

    def _checkPoint(self):

        if mpi.rank == 0:
            stepString = str(self.modelStep).zfill(5)
            stepOutputPath = self.outputPath + "/" + stepString
            newtime = sca.dimensionalise(self.modelTime, sca.units.year).magnitude
            try:
                os.mkdir(stepOutputPath)
            except FileExistsError:
                pass
            try:
                os.mkdir(stepOutputPath + "/h5")
            except FileExistsError:
                pass
            try:
                os.mkdir(stepOutputPath + "/xdmf")
            except FileExistsError:
                pass
            with open(stepOutputPath + "/time.json", "w") as f:
                json.dump({"time": f"{newtime:.3e} yrs"}, f)

        for name, subMesh in self._meshVarForSaving:
            obj = getattr(self, name)
            self.checkPointManager.saveField(
                name=name,
                field=obj,
                handle=self.meshHandle,
                step=self.modelStep,
                time=self.modelTime,
            )
        for name in self._swarmVarForSaving:
            obj = getattr(self, name)
            self.checkPointManager.saveField(
                name=name,
                field=obj,
                handle=self.swarmHandle,
                step=self.modelStep,
                time=self.modelTime,
            )
        self.figureManager.saveParticleViscosity(
            self.swarm, self.viscosityFn, self.modelStep
        )
        self.figureManager.saveStrainRate(
            strainRate2ndInvariant=self.strainRateField,
            mesh=self.mesh,
            step=self.modelStep,
        )
        self.figureManager.saveStressField(
            self.mesh, self.projectedStressField, self.modelStep
        )
        self.figureManager.saveTemperatureField(
            self.mesh, self.temperature, self.modelStep
        )
        self.figureManager.saveVelocity(
            self.velocityField, self.mesh, self.swarm, self.viscosityFn, self.modelStep
        )
        self.figureManager.saveMaterialVariable(
            self.materialVariable, self.swarm, self.modelStep
        )

    @abstractmethod
    def run(self):
        pass
