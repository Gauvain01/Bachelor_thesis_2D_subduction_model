import json
from abc import ABC
from typing import Dict, List, Tuple

import numpy as np
from black import out
from underworld import function as fn
from underworld import mpi
from underworld.swarm import Swarm
from underworld.systems import SwarmAdvector


class TracerParticle:
    def __init__(
        self,
        startingPosition: Tuple,
        name: str,
    ) -> None:

        self._startinPosition = startingPosition
        self.name = name

        if mpi.rank == 0:
            self.coords: Dict = {}
        mpi.barrier()
        self.saveInitialCoord()

    # called by the model
    def _buildSwarmAndAdvector(self, velocityField, mesh, outputPath):
        self._swarm = Swarm(mesh)
        self.outputPath = outputPath
        self._advector = SwarmAdvector(velocityField, self._swarm)
        self._addParticleToSwarm()

    @property
    def swarm(self):
        return self._swarm

    @property
    def startingPosition(self):
        return self._startinPosition

    def _addParticleToSwarm(self):
        coord_array = np.array(object=self.startingPosition, ndmin=2)
        self.swarm.add_particles_with_coordinates(coord_array)

    def saveInitialCoord(self):
        if mpi.rank == 0:
            coordArray = self.startingPosition
            coord = (coordArray[0], coordArray[1])
            item = {"time": 0, "coord": coord}
            self.coords[0] = item
        mpi.barrier()

    def saveCoord(self, step: int, time: int):

        if len(self.swarm.data) != 0:
            data = self.swarm.data
            mpi.comm.send(data, 0, 28)

        mpi.barrier()
        if mpi.rank == 0:
            coordArray = mpi.comm.recv(tag=28)[0]
            print(coordArray)
            coord = (coordArray[0], coordArray[1])
            item = {"time": time, "coord": coord}
            self.coords[step] = item
        mpi.barrier()

    def writeCoordsToJson(self):
        if mpi.rank == 0:
            path = f"{self.outputPath}/tracer_{self.name}_data.json"
            with open(path, "w") as f:
                json.dump(self.coords, f)
        mpi.barrier()

    @property
    def advector(self):
        return self._advector

    def advect(self, dt):
        self.advector.integrate(dt)
