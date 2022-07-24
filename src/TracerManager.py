from typing import List

from TracerParticle import TracerParticle


class TracerManager:
    def __init__(self) -> None:
        self._tracers: List[TracerParticle] = []

    @property
    def tracers(self):
        return self._tracers

    def addTracer(self, tracerParticle: TracerParticle):
        self.tracers.append(tracerParticle)

    def writeTracerData(self):
        for tracer in self.tracers:
            tracer.writeCoordsToJson()

    def saveCoordinatesForTracers(self, step, time):
        for tracer in self.tracers:
            tracer.saveCoord(step, time)

    def advectTracers(self, dt):
        for tracer in self.tracers:
            tracer.advect(dt)
