import UWGeodynamics as geo
from UWGeodynamics import visualisation as vis

u = geo.UnitRegistry

modelLength = 8000e3 * u.meter
modelDepth = 660e3 * u.meter


model = geo.Model(
    elementRes=(192, 48),
    minCoord=(0.0, 0.0),
    maxCoord=(4.0, 1.0),
    gravity=(0.0, -1.0),
    periodic=(True, False),
)

model.outputDir = "./src/output"


upperMantleShape = geo.shapes.Layer(top=model.top, bottom=0.4)
lowerMantleShape = geo.shapes.Layer(top=0.4, bottom=model.bottom)
lowerSlabShape = geo.shapes.Polygon(
    [
        (1.2, 0.925),
        (3.25, 0.925),
        (3.20, 0.900),
        (1.2, 0.900),
        (1.02, 0.825),
        (1.02, 0.850),
    ]
)
coreSlabShape = geo.shapes.Polygon(
    [
        (1.2, 0.975),
        (3.35, 0.975),
        (3.25, 0.925),
        (1.2, 0.925),
        (1.02, 0.850),
        (1.02, 0.900),
    ]
)
upperSlabShape = geo.shapes.Polygon(
    [
        (1.2, 1.000),
        (3.40, 1.000),
        (3.35, 0.975),
        (1.2, 0.975),
        (1.02, 0.900),
        (1.02, 0.925),
    ]
)


upperMantle = model.add_material(name="Upper Mantle", shape=upperMantleShape)
lowerMantle = model.add_material(name="Lower Mantle", shape=lowerMantleShape)
upperSlab = model.add_material(name="Upper Slab", shape=upperSlabShape)
lowerSlab = model.add_material(name="Lower Slab", shape=lowerSlabShape)
coreSlab = model.add_material(name="Core Slab", shape=coreSlabShape)


upperMantle.viscosity = 1.0
lowerMantle.viscosity = 100.0
upperSlab.viscosity = 500.0
lowerSlab.viscosity = 500.0
coreSlab.viscosity = 500.0


upperMantle.density = 0.0
lowerMantle.density = 0.0
upperSlab.density = 1.0
lowerSlab.density = 1.0
coreSlab.density = 1.0

upperSlab.plasticity = geo.VonMises(cohesion=0.06)
lowerSlab.plasticity = geo.VonMises(cohesion=0.06)


Fig = vis.Figure(figsize=(1200, 400))
Fig.Points(
    model.swarm, model.materialField, fn_size=2.0, colours="white green red purple blue"
)


model.set_velocityBCs(bottom=[0.0, 0.0], top=[None, 0.0])


model.run_for(nstep=2, checkpoint_interval=1)
Fig = vis.Figure(figsize=(1200, 400))
Fig.Points(
    model.swarm, model.materialField, fn_size=2.0, colours="white green red purple blue"
)
Fig.show()
