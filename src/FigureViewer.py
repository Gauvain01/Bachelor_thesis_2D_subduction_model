from underworld import visualisation

from FigureManager import FigureManager
from SubductionModel import SubductionModel


class FigureViewer:
    def __init__(self, modelPath) -> None:
        self.modelPath = modelPath

    def fromDb(self, step):
        dbPath = self.modelPath + "/FigStore.gldb"
        viewer = visualisation.Viewer(dbPath)

        viewer.step = step
        viewer.showall()

    def fromCheckPoint(self, step: int, model: SubductionModel):
        model.fromCheckpoint = True
        model._fromCheckpointStep = step
        model._initFromCheckPoint(step)
        mesh = model.mesh
        swarm = model.swarm
        figManager = FigureManager(model.outputPath, model.name, True)
        figManager.getParticlePlot(swarm, model.materialVariable)
        figManager.saveParticleViscosity(swarm, model.viscosityFn)
        figManager.saveStrainRate(
            model.rheologyCalculations.getStrainRateSecondInvariant(
                model.velocityField
            ),
            mesh,
        )
        figManager.saveTemperatureField(mesh, model.temperatureField)
        figManager.saveStress2ndInvariant(swarm, model.stress2ndInvariant)
        figManager.saveVelocity(model.velocityField, mesh, swarm, model.viscosityFn)
