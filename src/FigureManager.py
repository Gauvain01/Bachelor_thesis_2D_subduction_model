from underworld import function as fn
from underworld import mpi, visualisation

from CheckPointManager import CheckPointManager


class FigureManager:
    def __init__(self, outputPath, modelName, directView=False) -> None:
        self.name = modelName
        self.outputPath = outputPath
        self.store = visualisation.Store(f"{self.outputPath}/FigStore")
        self.directView = directView

    def saveFig(self, fig):
        if self.directView:
            fig.show()
        else:
            fig.save()

    def _getFig(self, title) -> visualisation.Figure:
        fig = visualisation.Figure(
            store=self.store,
            figsize=(1080, 720),
            title=title,
        )
        return fig

    def getParticlePlot(self, swarm, materialVariable):
        fig = self._getFig(f"{self.name} particles")
        fig.append(
            visualisation.objects.Points(
                swarm,
                materialVariable,
                pointSize=2,
                colours="green red purple blue yellow white orange",
                colourBar=False,
            )
        )
        self.saveFig(fig)

    def saveVelocity(self, velocityField, mesh, swarm, viscosityFn) -> None:
        fig = self._getFig(f"{self.name} Velocity")
        fig.append(
            visualisation.objects.Points(swarm, viscosityFn, pointSize=2, logScale=True)
        )
        fig.append(visualisation.objects.VectorArrows(mesh, velocityField))
        self.saveFig(fig)

    def saveStrainRate(self, strainRate2ndInvariant, mesh) -> None:
        fig = self._getFig(f"{self.name} Strain Rate 2nd Invariant")
        fig.append(
            visualisation.objects.Surface(
                mesh, strainRate2ndInvariant, onMesh=True, logScale=True
            )
        )
        self.saveFig(fig)

    def saveParticleViscosity(self, swarm, viscosityFn) -> None:
        fig = self._getFig(f"{self.name} Viscosity")

        fig.append(
            visualisation.objects.Points(swarm, viscosityFn, pointSize=2, logScale=True)
        )
        self.saveFig(fig)

    def saveTemperatureField(self, mesh, temperatureField):
        fig = self._getFig(f"{self.name} Viscosity")
        fig.append(visualisation.objects.Surface(mesh, temperatureField, onMesh=True))
        self.saveFig(fig)

    def saveStress2ndInvariant(self, swarm, stress2ndInvariant) -> None:
        fig = self._getFig(f"{self.name} Stress 2nd Invariant")

        fig.append(visualisation.objects.Points(swarm, stress2ndInvariant, pointSize=2))
        self.saveFig(fig)

    def incrementStoreStep(self) -> None:
        if mpi.rank == 0:
            self.store.step += 1
