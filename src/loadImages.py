from FigureViewer import FigureViewer

if __name__ == "__main__":

    def loadImages():
        testNumber = input("test number:")
        step = input("step:")
        lv = FigureViewer(f"src/output/test_{testNumber}")
        lv.fromDb(int(step))

    while True:
        loadImages()
        check = input("again y/n")
        if check == "n":
            break
