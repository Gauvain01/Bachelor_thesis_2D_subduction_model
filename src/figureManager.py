from underworld import function as fn
from underworld import visualisation


class FigureManager:
    def __init__(self, outputPath, modelName) -> None:
        self.name = modelName
        self.outputPath = outputPath
        self.store = visualisation.Store(f"{self.outputPath}/FigStore")

    def _getFig(self, title) -> visualisation.Figure:
        fig = visualisation.Figure(
            store=self.store,
            figsize=(960, 300),
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
        fig.save()

    def saveVelocityMagnitude(self, velocityField, mesh) -> None:
        fig = self._getFig(f"{self.name} Velocity Magnitude")
        fig.append(
            visualisation.objects.Surface(
                mesh,
                fn.math.sqrt(fn.math.dot(velocityField, velocityField)),
                onMesh=True,
            )
        )
        fig.save()

    def saveStrainRate(self, strainRate2ndInvariant, mesh) -> None:
        fig = self._getFig(f"{self.name} Strain Rate 2nd Invariant")
        fig.append(
            visualisation.objects.Surface(
                mesh, strainRate2ndInvariant, onMesh=True, logScale=True
            )
        )
        fig.save()

    def saveParticleViscosity(self, swarm, viscosityFn) -> None:
        fig = self._getFig(f"{self.name} Viscosity")

        fig.append(visualisation.objects.Points(swarm, viscosityFn, pointSize=2))
        fig.save()

    def saveStress2ndInvariant(self, swarm, stress2ndInvariant) -> None:
        fig = self._getFig(f"{self.name} Stress 2nd Invariant")

        fig.append(visualisation.objects.Points(swarm, stress2ndInvariant, pointSize=2))

    def incrementStoreStep(self) -> None:
        self.store.step += 1
