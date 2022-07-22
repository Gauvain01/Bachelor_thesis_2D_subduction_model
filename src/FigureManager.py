from underworld import function as fn
from underworld import mpi, visualisation

from CheckPointManager import CheckPointManager


class FigureManager:
    def __init__(self, outputPath, modelName, directView=False) -> None:
        self.name = modelName
        self.outputPath = outputPath
        self.store = visualisation.Store(f"{self.outputPath}/FigStore")
        self.directView = directView

    def saveFig(self, fig: visualisation.Figure):
        if self.directView:
            fig.show()
        else:
            fig.save()

    def _getFig(self, title, name) -> visualisation.Figure:
        fig = visualisation.Figure(
            store=self.store,
            figsize=(1080, 720),
            title=title,
            name=name + str(mpi.rank),
        )
        return fig

    def getParticlePlot(self, swarm, materialVariable, step):
        fig = self._getFig(f"{self.name} particles", f"particlePlot_{step}_")
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

    def saveVelocity(self, velocityField, mesh, swarm, materialVar, step) -> None:
        fig = self._getFig(f"{self.name} Velocity", f"velocity_{step}_")
        velocityMagnitude = fn.math.sqrt(fn.math.dot(velocityField, velocityField))
        fig.append(
            visualisation.objects.Points(
                swarm, velocityMagnitude, pointSize=2, logScale=True
            )
        )
        fig.append(
            visualisation.objects.VectorArrows(mesh, velocityField, autoscale=True)
        )
        fig.append(visualisation.objects.Points(swarm, materialVar))

        self.saveFig(fig)

    def saveStrainRate(self, strainRate2ndInvariant, mesh, step) -> None:
        fig = self._getFig(
            f"{self.name} Strain Rate 2nd Invariant", f"strainRate_{step}_"
        )
        fig.append(
            visualisation.objects.Surface(
                mesh, strainRate2ndInvariant, onMesh=True, logScale=True
            )
        )
        self.saveFig(fig)

    def saveMaterialVariable(self, materialVar, swarm, step) -> None:
        fig = self._getFig(f"{self.name} Material", f"materialVar_{step}_")
        fig.append(visualisation.objects.Points(swarm, materialVar))
        self.saveFig(fig)

    def saveParticleViscosity(self, swarm, viscosityFn, step) -> None:
        fig = self._getFig(f"{self.name} Viscosity", f"viscosity_{step}_")

        fig.append(
            visualisation.objects.Points(swarm, viscosityFn, pointSize=2, logScale=True)
        )
        self.saveFig(fig)

    def saveTemperatureField(self, mesh, temperatureField, step):
        fig = self._getFig(f"{self.name} Temperature", f"temperature_{step}_")
        fig.append(visualisation.objects.Surface(mesh, temperatureField))
        self.saveFig(fig)

    def saveStressField(self, mesh, stressField, step) -> None:
        fig = self._getFig(f"{self.name} StressField", f"stressField_{step}_")
        print(type(stressField))
        fig.append(visualisation.objects.Surface(mesh, stressField))
        self.saveFig(fig)

    def incrementStoreStep(self) -> None:
        if mpi.rank == 0:
            self.store.step += 1
