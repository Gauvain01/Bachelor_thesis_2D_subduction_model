from underworld import conditions
from underworld import function as fn
from underworld import mesh, mpi, swarm, visualisation

from modelParameters._Model_parameter_map import ModelParameterMap
from PlatePolygons import SubductionZonePolygons
from RheologyFunctions import RheologyFunctions


class SubductionModel:
    def __init__(
        self,
        name: str,
        modelParameterMap: ModelParameterMap,
        subductionZonePolygons: SubductionZonePolygons,
        stepAmount,
        stepTime,
        plateVelocity,
    ) -> None:
        self.name = name
        self.subductionZonePolygons = subductionZonePolygons
        self.parameters = modelParameterMap
        self.mesh = mesh.FeMesh_Cartesian(
            elementType="Q1/dQ0",
            elementRes=(512, 206),
            minCoord=(0.0, 0.0),
            maxCoord=(
                self.parameters.modelLength.nonDimensionalValue.magnitude,
                self.parameters.modelHeight.nonDimensionalValue.magnitude,
            ),
            periodic=[True, False],
        )
        self.rank = mpi.rank
        self.outputPath = "./output"

        self.swarm = swarm.Swarm(mesh=self.mesh)
        self.materialVariable = self.swarm.add_variable(dataType="int", count=1)
        self.previousStress = self.swarm.add_variable(dataType="double", count=3)
        self.dissipation = self.swarm.add_variable(dataType="double", count=1)
        self.storedEnergyRate = self.swarm.add_variable(dataType="double", count=1)
        self.swarmLayout = swarm.layouts.PerCellSpaceFillerLayout(
            swarm=self.swarm, particlesPerCell=20
        )
        self.swarm.populate_using_layout(self.swarmLayout)
        self.figStore = visualisation.Store("output/subduction")

        self.totalSteps = (stepAmount,)
        self.stepTime = self.parameters.scalingCoefficient.nonDimensionalizeUnderworld(
            stepTime.to_base_units()
        ).magnitude
        self.plateVelocity = (
            self.parameters.scalingCoefficient.nonDimensionalizeUnderworld(
                plateVelocity.to_base_units()
            ).magnitude
        )

        self._initializeFields()

        self.rheologyCalculations = RheologyFunctions(
            self.parameters, self.velocityField
        )

        self._initializeMaterialVariableData()
        self._initializePolygons()
        self._assignMaterialToParticles()
        self._initializeFields()
        self._setBoundaryConditions()
        self._initializeViscosityMapFn()
        self._initializeStressMapFn()

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
        self.lowerMantleIndex = 4
        self.coreSlabIndex = 3

    def _initializePolygons(self):
        upperSlabShape = self.subductionZonePolygons.getUpperSlabShapeArray()
        lowerSlabShape = self.subductionZonePolygons.getLowerSlabShapeArray()
        middleSlabShape = self.subductionZonePolygons.getMiddleSlabShapeArray()

        self.slabUpperPoly = fn.shape.Polygon(upperSlabShape)
        self.slabLowerPoly = fn.shape.Polygon(lowerSlabShape)
        self.slabCorePoly = fn.shape.Polygon(middleSlabShape)

    def _assignMaterialToParticles(self):
        self.materialVariable.data[:] = self.upperMantleIndex
        for index in range(len(self.swarm.particleCoordinates.data)):
            coord = self.swarm.particleCoordinates.data[index][:]
            if coord[1] < self.parameters.lowerMantleHeigth.nonDimensionalValue:
                self.materialVariable.data[index] = self.lowerMantleIndex
            if self.slabUpperPoly.evaluate(tuple(coord)):
                self.materialVariable.data[index] = self.upperSlabIndex
            if self.slabCorePoly.evaluate(tuple(coord)):
                self.materialVariable.data[index] = self.coreSlabIndex
            if self.slabLowerPoly.evaluate(tuple(coord)):
                self.materialVariable.data[index] = self.lowerSlabIndex

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
        visTopLayer = (
            self.rheologyCalculations.getEffectiveViscosityOfUpperLayerVonMises()
        )
        visCoreLayer = (
            self.rheologyCalculations.getEffectiveViscosityOfViscoElasticCore()
        )
        visBottomLayer = self.parameters.spBottomLayerViscosity.nonDimensionalValue

        viscosityMap = {
            self.upperMantleIndex: self.parameters.upperMantleViscosity.nonDimensionalValue,
            self.lowerMantleIndex: self.parameters.lowerMantleViscosity.nonDimensionalValue,
            self.lowerSlabIndex: visBottomLayer,
            self.coreSlabIndex: visCoreLayer,
            self.upperSlabIndex: fn.exception.SafeMaths(
                fn.misc.min(
                    visTopLayer,
                    self.parameters.spTopLayerViscosity.nonDimensionalValue,
                )
            ),
        }
        self.viscosityFn = fn.branching.map(
            fn_key=self.materialVariable, mapping=viscosityMap
        )

    def _initializeStressMapFn(self):
        Te = self.rheologyCalculations.getSymmetricStrainRateTensor()
        visFn = self.viscosityFn

        viscousStressFn = 2.0 * visFn * Te

        visEff = self.rheologyCalculations.getEffectiveViscosityOfViscoElasticCore()
        shearMod = self.parameters.coreShearModulus.nonDimensionalValue
        dt = self.parameters.timeScaleStress.nonDimensionalValue

        elasticStressFn = visEff / (shearMod * dt) * self.previousStress

        viscoElasticStressFn = viscousStressFn + elasticStressFn

        stressMap = {
            self.upperMantleIndex: viscousStressFn,
            self.lowerMantleIndex: viscousStressFn,
            self.lowerSlabIndex: viscousStressFn,
            self.coreSlabIndex: viscoElasticStressFn,
            self.upperSlabIndex: viscousStressFn,
        }
        self.stressFn = fn.branching.map(
            fn_key=self.materialVariable, mapping=stressMap
        )
        self.stress2ndInvariant = fn.tensor.second_invariant(self.stressFn)
