from __future__ import annotations

from abc import abstractmethod

from underworld import function as fn
from underworld import mesh, swarm, systems
from underworld.function import Function

from meshProjector import meshProjector
from modelParameters import ModelParameterMap


class BaseModel:
    def __init__(
        self, modelParameters: ModelParameterMap, totalSteps: int, checkPointSteps: int
    ) -> None:
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
        self.currentStep = 0
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

    def _addSwarmVariable(self, name: str, dataType: str, count: int):
        self.swarm: swarm.Swarm
        field = self.swarm.add_variable(dataType, count)
        setattr(self, name, field)

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

    @abstractmethod
    def _checkPoint(self):
        pass

    @abstractmethod
    def run(self):
        pass
