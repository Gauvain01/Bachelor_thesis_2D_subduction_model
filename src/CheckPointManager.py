import json
import os
from typing import Union

from underworld import mesh as Mesh
from underworld import mpi, swarm
from underworld.swarm import Swarm


class CheckPointManager:
    def __init__(self, modelName, outputPath) -> None:
        self.ModelName = modelName
        self.outputPath = outputPath

    def _getStepOutputPath(self, step: int):
        stepString = str(step).zfill(5)
        stepOutputPath = self.outputPath + "/" + stepString
        return stepOutputPath

    def loadField(
        self, name: str, step: int, field: Union[Mesh.MeshVariable, swarm.SwarmVariable]
    ):
        with mpi.call_pattern("collective"):
            path = self._getH5Path(step) + name + ".h5"
            field.load(path, interpolate=True)

    def saveField(
        self,
        name: str,
        field: Union[Mesh.MeshVariable, swarm.SwarmVariable],
        handle,
        step: int,
        time,
    ):
        with mpi.call_pattern("collective"):
            h5Path = self._getH5Path(step)
            xdmfPath = self._getXdmfPath(step)

            fieldHnd = field.save(h5Path + f"{name}.h5")
            if isinstance(field, Mesh.MeshVariable):
                field.xdmf(
                    filename=xdmfPath + f"{name}.xdmf",
                    fieldSavedData=fieldHnd,
                    varname=name,
                    meshSavedData=handle,
                    meshname="mesh",
                    modeltime=time,
                )
            if isinstance(field, swarm.SwarmVariable):
                field.xdmf(
                    filename=xdmfPath + f"{name}.xdmf",
                    varSavedData=fieldHnd,
                    varname=name,
                    swarmSavedData=handle,
                    swarmname="swarm",
                    modeltime=time,
                )

    def _getH5Path(self, step: int):
        return self._getStepOutputPath(step) + "/h5/"

    def _getXdmfPath(self, step: int):
        return self._getStepOutputPath(step) + "/xdmf/"

    def getMesh(self, mesh) -> Mesh.FeMesh_Cartesian:
        meshPath = self.outputPath + "mesh.00000.h5"
        mesh.load(filename=meshPath)

    def getSwarm(self, swarm, step) -> Swarm:
        swarmPath = self._getH5Path(step) + "swarm.h5"
        swarm.load(swarmPath)

    def getMaterialVariable(self, step, swarm: Swarm):
        mvarPath = self._getH5Path(step) + "materialVariable.h5"
        materialVar = swarm.add_variable(dataType="int", count=1)
        materialVar.load(mvarPath)
        return materialVar

    def getPreviousStress(self, step, swarm: Swarm):
        path = self._getH5Path(step) + "previousStress.h5"
        var = swarm.add_variable(dataType="double", count=3)
        var.load(path)
        return var

    def getVelocityField(self, step, mesh: Mesh.FeMesh_Cartesian):
        path = self._getH5Path(step) + "velocityField.h5"
        field = Mesh.MeshVariable(mesh=mesh, nodeDofCount=2)
        field.load(path, interpolate=True)
        return field

    def getPressureField(self, step, mesh: Mesh.FeMesh_Cartesian):
        path = self._getH5Path(step) + "pressureField.h5"
        field = Mesh.MeshVariable(mesh=mesh, nodeDofCount=1)
        field.load(path, interpolate=True)
        return field

    def getTemperatureField(self, step, mesh: Mesh.FeMesh_Cartesian):
        path = self._getH5Path(step) + "temperatureField.h5"
        field = Mesh.MeshVariable(mesh=mesh, nodeDofCount=1)
        field.load(path, interpolate=True)
        return field

    def getTemperatureDotField(self, step, mesh: Mesh.FeMesh_Cartesian):
        path = self._getH5Path(step) + "temperatureDotField.h5"
        field = Mesh.MeshVariable(mesh=mesh, nodeDofCount=1)
        field.load(path)
        return field

    def getLastTime(self, step):
        path = self._getStepOutputPath(step) + "/time.json"
        with open(path, "r") as f:
            return json.load(f)

    def checkPoint(
        self,
        *,
        step,
        swarm: Swarm,
        mesh,
        materialVariable: Mesh.MeshVariable,
        previousStress: swarm.SwarmVariable,
        velocityField,
        pressureField,
        temperatureField,
        temperatureDotField,
        figureManager,
        meshHandle,
        time,
        strainRate2ndInvariant,
        viscosityFn,
        stress2ndInvariant,
    ):
        stepString = str(step).zfill(5)
        stepOutputPath = self.outputPath + "/" + stepString

        if mpi.rank == 0:
            os.mkdir(stepOutputPath)
            os.mkdir(stepOutputPath + "/h5")
            os.mkdir(stepOutputPath + "/xdmf")

        mpi.barrier()

        h5Path = stepOutputPath + "/h5/"
        xdmfPath = stepOutputPath + "/xdmf/"

        swarmHnd = swarm.save(h5Path + "swarm.h5")
        materialVariableHnd = materialVariable.save(h5Path + "materialVariable.h5")
        previousStressHnd = previousStress.save(h5Path + "previousStress.h5")
        temperatureDotHnd = temperatureDotField.save(
            h5Path + "temperatureDotField" + ".h5", meshHandle
        )
        velocityHnd = velocityField.save(h5Path + "velocityField" + ".h5", meshHandle)
        pressureHnd = pressureField.save(h5Path + "pressureField" + ".h5", meshHandle)
        temperatureHnd = temperatureField.save(
            h5Path + "temperatureField" + ".h5", meshHandle
        )

        materialVariable.xdmf(
            xdmfPath + "materialVariable.xdmf",
            varSavedData=materialVariableHnd,
            varname="materialVariable",
            swarmSavedData=swarmHnd,
            swarmname="swarm",
            modeltime=time,
        )
        previousStress.xdmf(
            filename=xdmfPath + "previousStress.xdmf",
            varSavedData=previousStressHnd,
            varname="previousStress",
            swarmSavedData=swarmHnd,
            swarmname="swarm",
            modeltime=time,
        )
        velocityField.xdmf(
            filename=xdmfPath + "velocityField.xdmf",
            fieldSavedData=velocityHnd,
            varname="velocityField",
            meshname="mesh",
            meshSavedData=meshHandle,
            modeltime=time,
        )
        temperatureField.xdmf(
            filename=xdmfPath + "temperatureField.xdmf",
            fieldSavedData=temperatureHnd,
            varname="temperatureField",
            meshSavedData=meshHandle,
            meshname="mesh",
            modeltime=time,
        )
        temperatureDotField.xdmf(
            filename=xdmfPath + "temperatureDotField.xdmf",
            fieldSavedData=temperatureDotHnd,
            varname="temperatureDotField",
            meshSavedData=meshHandle,
            meshname="mesh",
            modeltime=time,
        )
        pressureField.xdmf(
            filename=xdmfPath + "pressureField.xdmf",
            fieldSavedData=pressureHnd,
            varname="pressureField",
            meshSavedData=meshHandle,
            meshname="mesh",
            modeltime=time,
        )

        # figureManager.saveParticleViscosity(swarm, viscosityFn)
        # figureManager.saveStrainRate(strainRate2ndInvariant, mesh)
        # figureManager.saveStress2ndInvariant(swarm, stress2ndInvariant)
        # figureManager.saveVelocity(velocityField, mesh, swarm, viscosityFn)
