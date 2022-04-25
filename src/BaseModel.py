from __future__ import annotations

from abc import abstractmethod
from typing import Union

from underworld import function as fn
from underworld import mesh, mpi, swarm, systems
from underworld.function import Function

from CheckPointManager import CheckPointManager
from meshProjector import meshProjector
from modelParameters import ModelParameterMap


class BaseModel:
    def __init__(
        self,
        modelParameters: ModelParameterMap,
        totalSteps: int,
        checkPointSteps: int,
        name: int,
    ) -> None:
        self.name = name
        self.totalSteps = totalSteps
        self.checkPointSteps = checkPointSteps
        self.parameters = modelParameters
        self._setMesh()
        self._setSwarm()
        self._addMeshVariable("pressureField", "double", 1, subMesh=True)
        self._addMeshVariable("velocityField", "double", 2)
        self._addMeshVariable("_strainRateField", "double", 1, subMesh=True)
        self._addSwarmVariable("_viscosityField", "double", 1)
        self._addMeshVariable("_projectedViscosity", "double", 1, subMesh=True)
        self._addSwarmVariable("_stressField", "double", 1)
        self._addMeshVariable(
            "_projectedStress", "double", nodeDofCount=1, subMesh=True
        )
        self._addSwarmVariable("_stressTensor", "double", count=3)
        self._addMeshVariable(
            "_projectedStressTensor", "double", nodeDofCount=3, subMesh=True
        )
        self._addSwarmVariable("_materialVariable", dataType="int", count=1)
        self.modelStep = 0
        self.modelTime = 0
        self._solutionExists = fn.misc.constant(False)
        self._hasTemperatureDotField = False
        self._hasTemperatureField = False

    def _addMeshVariable(
        self, name: str, dataType: str, nodeDofCount, subMesh: bool = False
    ):
        if subMesh:
            field = mesh.MeshVariable(self.mesh.subMesh, nodeDofCount, dataType)
        else:
            field = mesh.MeshVariable(self.mesh, nodeDofCount, dataType)
        setattr(self, name.value, field)
        self._meshVarForSaving.append(name)

    def _addSwarmVariable(self, name: str, dataType: str, count: int):
        self.swarm: swarm.Swarm
        field = self.swarm.add_variable(dataType, count)
        setattr(self, name, field)
        self._swarmVarForSaving.append(name)

    @abstractmethod
    def _init_temperature_variables(self):
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

    @property
    @abstractmethod
    def materialVariable(self):
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
    def _swarmVarForSaving(self):
        obj = []
        return obj

    @property
    def _meshVarForSaving(self):
        obj = []
        return obj

    @property
    def checkPointManager(self) -> CheckPointManager:
        obj = CheckPointManager(self.name, self.outputPath)
        return obj

    @property
    def temperature(self):
        if not self._hasTemperature:
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
        stress = fn.tensor.second_invariant(self.stressFn)
        self._stressField.data[...] = stress.evaluate(self.swarm)
        meshProjector(self._projectedStress, self._stressField)
        return self._projectedStressField

    @property
    def projectedStressTensor(self):
        self._stressTensor.data[...] = self.stressFn.evaluate(self.swarm)
        meshProjector(self._projectedStressTensor, self._stressTensor)
        return self._projectedStressTensor

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
        return self._stokes

    @property
    def swarmAdvector(self):
        obj = systems.SwarmAdvector(self.velocityField, self.swarm, order=2)
        return obj

    @property
    def solver(self):
        self._solver = systems.Solver(self.stokes)
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

    def solve(self):
        self.solver.solve(
            nonLinearIterate=True, nonLinearTolerance=1e-2, print_stats=True
        )

    @abstractmethod
    def _update(self):
        pass

    def _checkPoint(self):
        for name in self._meshVarForSaving:
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

    @abstractmethod
    def _run(self):
        pass

    def _cleanup(self):
        pass

    def run(self):
        self._run()
        self._cleanup()
