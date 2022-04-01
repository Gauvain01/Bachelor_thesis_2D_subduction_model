import numpy as np
import UWGeodynamics as geo
from UWGeodynamics import visualisation as vis

from model_parameter_sets.Strak_2021_model_parameters import strak2021ModelParameterSet
from PlatePolygons import OceanicPLatePolygons

parameterSet = strak2021ModelParameterSet
# nonDimensionalizeParameters
parameterSet.nonDimensionalizeParameters()
import math

u = geo.UnitRegistry


model = geo.Model(
    elementRes=(1024, 512),
    minCoord=(0.0, 0.0),
    maxCoord=(
        parameterSet.modelLength.nonDimensionalValue.magnitude,
        parameterSet.modelHeight.nonDimensionalValue.magnitude,
    ),
    gravity=(0.0, -1.0),
    periodic=(False, False),
)

pp = OceanicPLatePolygons(
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
model.outputDir = "./src/output"


upperMantleShape = geo.shapes.Layer(top=model.top, bottom=model.bottom)

lowerSlabShape = geo.shapes.Polygon(pp.getLowerSlabShapeArray())
coreSlabShape = geo.shapes.Polygon(pp.getMiddleSlabShapeArray())
upperSlabShape = geo.shapes.Polygon(pp.getUpperSlabShapeArray())
crustShape = geo.shapes.Polygon(pp.getCrustSlabShapeArray())
foreArcShape = geo.shapes.Polygon(pp.getLithosphericMantleShapeForeArc())
backArcShape = geo.shapes.Polygon(pp.getLithosphericMantleShapeFarBackArc())

upperMantle = model.add_material(name="Upper Mantle", shape=upperMantleShape)

upperSlab = model.add_material(name="Upper Slab", shape=upperSlabShape)
lowerSlab = model.add_material(name="Lower Slab", shape=lowerSlabShape)
coreSlab = model.add_material(name="Core Slab", shape=coreSlabShape)
crust = model.add_material(name="crust", shape=crustShape)
foreArc = model.add_material(name="foreArc", shape=foreArcShape)
backArc = model.add_material(name="backArc", shape=backArcShape)


crust.viscosity = 1.0
foreArc.viscosity = 1.0
backArc.viscosity = 1.0
upperMantle.viscosity = 1.0
upperSlab.viscosity = 500.0
lowerSlab.viscosity = 500.0
coreSlab.viscosity = 500.0

crust.density = 0.0
foreArc.density = 0.0
backArc.density = 0.0
upperMantle.density = 0.0
upperSlab.density = 1.0
lowerSlab.density = 1.0
coreSlab.density = 1.0

# upperSlab.plasticity = geo.VonMises(cohesion=0.06)
# lowerSlab.plasticity = geo.VonMises(cohesion=0.06)


Fig = vis.Figure(figsize=(2160, 1080))
Fig.Points(
    model.swarm,
    model.materialField,
    fn_size=2.0,
    colours="green red purple blue yellow white orange",
)
Fig.save("figure_1.png")
Fig.show()

# model.set_velocityBCs(bottom=[0.0, 0.0], top=[None, 0.0])


# model.run_for(nstep=2, checkpoint_interval=1)
# Fig = vis.Figure(figsize=(1200, 400))
# Fig.Points(
#     model.swarm, model.materialField, fn_size=2.0, colours="white green red purple blue"
# )
# Fig.show()
