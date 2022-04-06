from __future__ import annotations

import os

from underworld import conditions
from underworld import function as fn
from underworld import mesh, mpi, swarm, systems
from underworld import underworld as uw
from underworld import visualisation

from model_parameter_sets.Strak_2021_model_parameters import strak2021ModelParameterSet
from model_parameters.Model_parameter_set import ModelParameterSet
from PlatePolygons import SubductionZonePolygons
from rheologyAlgorithms import RheologyCalculations

u = uw.scaling.units
# model = geo.Model(
#     elementRes=(1024, 206),
#     minCoord=(0.0, 0.0),
#     maxCoord=(
#         parameterSet.modelLength.nonDimensionalValue.magnitude,
#         parameterSet.modelHeight.nonDimensionalValue.magnitude,
#     ),
#     periodic=(True, False),
# )


# model.outputDir = "./src/output"

# upperMantleShape = geo.shapes.Layer(top=model.top, bottom=model.bottom)
# lowerSlabPP = pp.getLowerSlabShapeArray()
# coreSlabPP = pp.getMiddleSlabShapeArray()
# upperSlabPP = pp.getUpperSlabShapeArray()
# crustSlabPP = pp.getCrustSlabShapeArray()
# foreArcPP = pp.getLithosphericMantleShapeForeArc()
# backArcPP = pp.getLithosphericMantleShapeFarBackArc()

# lowerSlabShape = geo.shapes.Polygon(lowerSlabPP)
# coreSlabShape = geo.shapes.Polygon(coreSlabPP)
# upperSlabShape = geo.shapes.Polygon(upperSlabPP)
# crustShape = geo.shapes.Polygon(crustSlabPP)
# foreArcShape = geo.shapes.Polygon(foreArcPP)
# backArcShape = geo.shapes.Polygon(backArcPP)

# upperMantle = model.add_material(name="Upper Mantle", shape=upperMantleShape)

# upperSlab = model.add_material(name="Upper Slab", shape=upperSlabShape)
# lowerSlab = model.add_material(name="Lower Slab", shape=lowerSlabShape)
# coreSlab = model.add_material(name="Core Slab", shape=coreSlabShape)
# crust = model.add_material(name="crust", shape=crustShape)
# foreArc = model.add_material(name="foreArc", shape=foreArcShape)
# backArc = model.add_material(name="backArc", shape=backArcShape)


# crust.viscosity = 1.0
# foreArc.viscosity = 1.0
# backArc.viscosity = 1.0
# upperMantle.viscosity = 1.0
# upperSlab.viscosity = 500.0
# lowerSlab.viscosity = 500.0
# coreSlab.viscosity = 500.0

# crust.density = 0.0
# foreArc.density = 0.0
# backArc.density = 0.0
# upperMantle.density = 0.0
# upperSlab.density = 1.0
# lowerSlab.density = 1.0
# coreSlab.density = 1.0

# # upperSlab.plasticity = geo.VonMises(cohesion=0.06)
# # lowerSlab.plasticity = geo.VonMises(cohesion=0.06)


# Fig = vis.Figure(figsize=(2160, 1080), axis=True)
# Fig.Points(
#     model.swarm,
#     model.materialField,
#     fn_size=2.0,
#     colours="green red purple blue yellow white orange",
# )
# Fig.save("figure_2.png")
# Fig.show()

# # model.set_velocityBCs(bottom=[0.0, 0.0], top=[None, 0.0])


# # model.run_for(nstep=2, checkpoint_interval=1)
# # Fig = vis.Figure(figsize=(1200, 400))
# # Fig.Points(
# #     model.swarm, model.materialField, fn_size=2.0, colours="white green red purple blue"
# # )
# # Fig.show()


class SubductionModel:
    def __init__(
        self,
        name: str,
        modelParameterSet: ModelParameterSet,
        subductionZonePolygons: SubductionZonePolygons,
        stepAmount,
        stepTime,
        plateVelocity,
    ) -> None:
        self.name = name
        self.subductionZonePolygons = subductionZonePolygons
        self.parameterSet = modelParameterSet
        self.mesh = mesh.FeMesh_Cartesian(
            elementType="Q1/dQ0",
            elementRes=(1024, 512),
            minCoord=(0.0, 0.0),
            maxCoord=(
                self.parameterSet.modelLength.nonDimensionalValue.magnitude,
                self.parameterSet.modelHeight.nonDimensionalValue.magnitude,
            ),
            periodic=[True, False],
        )
        self.rank = mpi.rank
        self.outputPath = "./output"

        self.swarm = swarm.Swarm(mesh=self.mesh)
        self.materialVariable = self.swarm.add_variable(dataType="int", count=1)
        self.swarmLayout = swarm.layouts.PerCellSpaceFillerLayout(
            swarm=self.swarm, particlesPerCell=20
        )
        self.swarm.populate_using_layout(self.swarmLayout)
        self.figStore = visualisation.Store("output/subduction")

        self.totalSteps = (stepAmount,)
        self.stepTime = (
            self.parameterSet._scalingCoefficient.nonDimensionalizeUnderworld(
                stepTime.to_base_units()
            ).magnitude
        )
        self.plateVelocity = (
            self.parameterSet._scalingCoefficient.nonDimensionalizeUnderworld(
                plateVelocity.to_base_units()
            ).magnitude
        )

        self._initializeFields()

        self.rheologyCalculations = RheologyCalculations(
            self.parameterSet, self.velocityField
        )

        self._initializeMaterialVariableData()
        self._initializePolygons()
        self._assignMaterialToParticles()
        self._initializeFields()
        self._setBoundaryConditions()
        self._initializeViscosityMapFn()

    def _setBoundaryConditions(self):
        self.verticalWalls = (
            self.mesh.specialSets["Left_VertexSet"]
            + self.mesh.specialSets["Right_VertexSet"]
        )
        self.horizontalWalls = (
            self.mesh.specialSets["Top_VertexSet"]
            + self.mesh.specialSets["Bottom_VertexSet"]
        )
        self.rightWall = self.mesh.specialSets["Right_VertexSet"]
        self.leftWall = self.mesh.specialSets["Left_VertexSet"]
        self.topWall = self.mesh.specialSets["Top_VertexSet"]
        self.bottomWall = self.mesh.specialSets["Bottom_VertexSet"]

        fixed = self.mesh.specialSets["Empty"]
        fixNode = self.rightWall + self.bottomWall
        fixed += fixNode

        self.VelocityBoundaryCondition = conditions.DirichletCondition(
            variable=self.velocityField,
            indexSetsPerDof=(self.verticalWalls + self.horizontalWalls),
        )
        self.temperatureBoundaryCondition = conditions.DirichletCondition(
            variable=self.temperatureField,
            indexSetsPerDof=(self.topWall + self.verticalWalls),
        )

    def _initializeFields(self):
        self.velocityField = mesh.MeshVariable(mesh=self.mesh, nodeDofCount=2)
        self.pressureField = mesh.MeshVariable(mesh=self.mesh, nodeDofCount=1)
        self.temperatureField = mesh.MeshVariable(mesh=self.mesh, nodeDofCount=1)
        self.temperatureDotField = mesh.MeshVariable(mesh=self.mesh, nodeDofCount=1)
        self.pressureField.data[:] = 0.0
        self.temperatureDotField.data[:] = 0.0
        self.temperatureField.data[:] = 0

    def _initializeMaterialVariableData(self):
        self.upperMantleIndex = 0
        self.upperSlabIndex = 1
        self.lowerSlabIndex = 2
        self.crustIndex = 4
        self.lithoSphericMantleForeArcIndex = 5
        self.lithoSphericMantleBackArcIndex = 6
        self.middleSlabIndex = 3

    def _initializePolygons(self):
        upperSlabShape = self.subductionZonePolygons.getUpperSlabShapeArray()
        lowerSlabShape = self.subductionZonePolygons.getLowerSlabShapeArray()
        middleSlabShape = self.subductionZonePolygons.getMiddleSlabShapeArray()

        crustShape = self.subductionZonePolygons.getCrustSlabShapeArray()
        lithosphericMantleForeArcShape = (
            self.subductionZonePolygons.getLithosphericMantleShapeForeArc()
        )
        lithosphericMantleBackArcShape = (
            self.subductionZonePolygons.getLithosphericMantleShapeFarBackArc()
        )

        self.slabUpperPoly = fn.shape.Polygon(upperSlabShape)
        self.slabLowerPoly = fn.shape.Polygon(lowerSlabShape)
        self.slabMiddlePoly = fn.shape.Polygon(middleSlabShape)

        self.crustPoly = fn.shape.Polygon(crustShape)
        self.lithosphericMantleForeArcPoly = fn.shape.Polygon(
            lithosphericMantleForeArcShape
        )
        self.lithosphericMantleBackArcPoly = fn.shape.Polygon(
            lithosphericMantleBackArcShape
        )

    def _assignMaterialToParticles(self):
        self.materialVariable.data[:] = self.upperMantleIndex
        for index in range(len(self.swarm.particleCoordinates.data)):
            coord = self.swarm.particleCoordinates.data[index][:]
            if self.slabUpperPoly.evaluate(tuple(coord)):
                self.materialVariable.data[index] = self.upperSlabIndex
            if self.slabMiddlePoly.evaluate(tuple(coord)):
                self.materialVariable.data[index] = self.middleSlabIndex
            if self.slabLowerPoly.evaluate(tuple(coord)):
                self.materialVariable.data[index] = self.lowerSlabIndex

            if self.crustPoly.evaluate(tuple(coord)):
                self.materialVariable.data[index] = self.crustIndex
            if self.lithosphericMantleForeArcPoly.evaluate(tuple(coord)):
                self.materialVariable.data[index] = self.lithoSphericMantleForeArcIndex
            if self.lithosphericMantleBackArcPoly.evaluate(tuple(coord)):
                self.materialVariable.data[index] = self.lithoSphericMantleForeArcIndex

    def getParticlePlot(self):
        fig = visualisation.Figure(
            self.figStore, figsize=(2160, 1080), name=f"{self.name} particles"
        )
        fig.append(
            visualisation.objects.Points(
                self.swarm,
                self.materialVariable,
                pointSize=2,
                colours="green red purple blue yellow white orange",
                colourBar=False,
            )
        )
        fig.show()

    def _initializeViscosityMapFn(self):
        def _safeViscosity(func, visMin, visMax):
            return fn.misc.max(visMin, fn.misc.min(visMax, func))

        upperMantleYieldVis = _safeViscosity(
            self.rheologyCalculations.getEffectiveViscosityDiffusionCreep(),
            self.parameterSet.lowerViscosityCutOff.nonDimensionalValue,
            self.parameterSet.upperViscosityCutOff.nonDimensionalValue,
        )
        SPplateUpperYieldVis = fn.exception.SafeMaths(
            fn.misc.min(
                self.rheologyCalculations.getEffectiveViscosityOfUpperLayerVonMises(),
                (self.parameterSet.referenceViscosity.nonDimensionalValue * 0.1),
            )
        )

        SPplateCoreYieldVis = self.parameterSet.spCoreLayerViscosity.nonDimensionalValue
        SPplateLowerYieldVis = (
            self.parameterSet.spBottomLayerViscosity.nonDimensionalValue
        )

        OPcrustVis = self.parameterSet.opCrustLayerViscosity.nonDimensionalValue
        OPLithoForeBackVis = (
            self.parameterSet.opLithosphericMantleViscosityInForearcAndBackarc.nonDimensionalValue
        )
        OPlithoFarBackVis = (
            self.parameterSet.opLithosphericMantleViscosityInFarBackarc.nonDimensionalValue
        )

        viscosityMap = {
            self.upperMantleIndex: upperMantleYieldVis,
            self.upperSlabIndex: SPplateUpperYieldVis,
            self.middleSlabIndex: SPplateCoreYieldVis,
            self.crustIndex: OPcrustVis,
            self.lithoSphericMantleForeArcIndex: OPLithoForeBackVis,
            self.lithoSphericMantleBackArcIndex: OPlithoFarBackVis,
            self.lowerSlabIndex: SPplateLowerYieldVis,
        }

        self.viscosityMapFn = fn.branching.map(
            fn_key=self.materialVariable, mapping=viscosityMap
        )


parameterSet = strak2021ModelParameterSet
# nonDimensionalizeParameters
parameterSet.nonDimensionalizeParameters()

pp = SubductionZonePolygons(
    parameterSet,
    29,
    200e3 * u.meter,
    6000e3 * u.meter,
    30e3 * u.meter,
    20e3 * u.meter,
    30e3 * u.meter,
    100e3 * u.meter,
    500e3 * u.meter,
    200e3 * u.meter,
    3000e3 * u.meter,
    30e3 * u.meter,
    30e3 * u.meter,
    120e3 * u.meter,
)
model = SubductionModel("test", parameterSet, pp)

if __name__ == "__main__":
    outputPath = os.path.join(os.path.abspath("."), "output/")

    if uw.mpi.rank == 0:
        if not os.path.exists(outputPath):
            os.makedirs(outputPath)
    mpi.barrier()
    model.getParticlePlot()
