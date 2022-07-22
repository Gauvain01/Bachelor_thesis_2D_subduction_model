import os

from FigureViewer import FigureViewer

if __name__ == "__main__":

    def loadImages(testNumber: int, step):

        lv = FigureViewer(f"src/output/test_{testNumber}")
        lv.fromDb(int(step))

    def getStepNumbers():
        stepList = []
        while True:

            item = input("step use * to stop: ")
            if item == "*":
                break
            stepList.append(item)
        return stepList

    def main():

        inputNumber = input("test number: ")
        stepList = getStepNumbers()
        for number in stepList:
            loadImages(inputNumber, number)

        print(f"finished images for steps {stepList}")

    main()
