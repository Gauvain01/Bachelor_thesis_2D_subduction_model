from underworld import visualisation


class FigureViewer:
    def __init__(self, modelPath) -> None:
        self.modelPath = modelPath

    def fromDb(self, step: int):
        dbPath = self.modelPath + "/FigStore.gldb"
        viewer = visualisation.Viewer(dbPath)

        viewer.step = step
        for name in viewer.figures:
            name: str
            field, rank, stepValue = name.split("_")
            if rank == "0":
                viewer.figure(name)
                viewer.show(
                    filename=self.modelPath + f"/{str(step).zfill(5)}/{name}.png"
                )
