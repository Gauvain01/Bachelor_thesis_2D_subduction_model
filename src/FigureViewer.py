from underworld import visualisation


class FigureViewer:
    def __init__(self, modelPath) -> None:
        self.modelPath = modelPath

    def fromDb(self, step):
        dbPath = self.modelPath + "/FigStore.gldb"
        viewer = visualisation.Viewer(dbPath)

        viewer.step = step
        viewer.showall()
