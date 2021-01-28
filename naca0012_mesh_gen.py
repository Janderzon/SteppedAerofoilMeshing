import naca_4_series_points
import gmsh
import sys
import math
import aerofoil_scale

################################################################################
#   Input Parameters
################################################################################

n = 50                          #Number points on the aerofoil surface
surfaceOrder = 2
farFieldSize = 10               #Width and height of the far field volume boundary
farFieldMeshSize = 1            #Target mesh size at the far field volume boundary
boundaryLayerThickness = 0.2    #Thickness of the boundary layer
boundaryLayerNumCells = 20      #Number of cells in the boundary layer in the direction normal to the aerofoil surface
boundaryLayerProgression = 1.3  #Growth rate of the cells in the boundary layer
wakeLength = 10
wakeEndThickness = 0.2
wakeLengthPoints = 100

################################################################################
#   Aerofoil Geometry
################################################################################

#Get the coordinates of the aerofoil
coords = naca_4_series_points.points(0, 0, 12, n*surfaceOrder, "cos")

#Initialize Gmsh and name the model
gmsh.initialize()
gmsh.model.add("NACA0012")

#Create the aerofoil points
aerofoilPoints = []
for i in range(0, len(coords)):
    aerofoilPoints.append(gmsh.model.geo.addPoint(coords[i][0], coords[i][1], 0))

#Create the aerofoil lines
aerofoilCurves = []
for i in range(0, n*2):
    j = i*surfaceOrder
    curvePoints = []
    for k in range(0, surfaceOrder+1):
        curvePoints.append(aerofoilPoints[k+j])
    aerofoilCurves.append(gmsh.model.geo.addSpline(curvePoints))
    gmsh.model.geo.mesh.setTransfiniteCurve(aerofoilCurves[i], 2)
trailingEdgeLine = gmsh.model.geo.addLine(aerofoilPoints[0], aerofoilPoints[-1])
gmsh.model.geo.mesh.setTransfiniteCurve(trailingEdgeLine, 2)

#gmsh.model.geo.mesh.setTransfiniteCurve(aerofoilSpline, n)
#upperTrailingEdgeLine = gmsh.model.geo.addLine(trailingEdgePoint, aerofoilPoints[0])
#gmsh.model.geo.mesh.setTransfiniteCurve(upperTrailingEdgeLine, 3)
#lowerTrailingEdgeLine = gmsh.model.geo.addLine(aerofoilPoints[-1], trailingEdgePoint)
#gmsh.model.geo.mesh.setTransfiniteCurve(lowerTrailingEdgeLine, 3)

#Create the aerofoil curve loop
#aerofoilLoop = gmsh.model.geo.addCurveLoop([upperTrailingEdgeLine, aerofoilSpline, lowerTrailingEdgeLine])

################################################################################
#   Boundary Layer Geometry
################################################################################

#Create boundary layer points
boundaryLayerPoints = []
numStraightPoints = round((n+1)*0.8)
numCircPoints = (n-round((n+1)*0.8))*2
for i in range(0, numStraightPoints):
    j = i*surfaceOrder
    boundaryLayerPoints.append(gmsh.model.geo.addPoint(coords[j][0], boundaryLayerThickness, 0))
for i in range(0, numCircPoints+1):
    angle = i*math.pi/numCircPoints
    boundaryLayerPoints.append(gmsh.model.geo.addPoint(0.05-boundaryLayerThickness*math.sin(angle), boundaryLayerThickness*math.cos(angle), 0))
for i in range(numCircPoints+numStraightPoints+1, n*2+1):
    j = i*surfaceOrder
    boundaryLayerPoints.append(gmsh.model.geo.addPoint(coords[j][0], -boundaryLayerThickness, 0))

#Create the boundary layer lines
boundaryLayerLines = []
for i in range(0, len(boundaryLayerPoints)-1):
    boundaryLayerLines.append(gmsh.model.geo.addLine(boundaryLayerPoints[i], boundaryLayerPoints[i+1]))
    gmsh.model.geo.mesh.setTransfiniteCurve(boundaryLayerLines[i], 2)

#Create the boundary layer curve loop
# boundaryLayerLoop = gmsh.model.geo.addCurveLoop([upperRearBoundaryLayerLine, upperBoundaryLayerLine, leadingEdgeCircLine, lowerBoundaryLayerLine, lowerRearBoundaryLayerLine])
#gmsh.model.geo.mesh.setTransfiniteCurve(boundaryLayerLoop, n)

#Create boundary layer surface
#boundaryLayerSurface = gmsh.model.geo.addPlaneSurface([boundaryLayerLoop, aerofoilLoop])
#gmsh.model.geo.mesh.setTransfiniteSurface(boundaryLayerSurface)
#gmsh.model.geo.mesh.setRecombine(2,boundaryLayerSurface)

#Create lines to divide up boundary layer into quadrilaterals
boundaryLayerDividerLines = []
for i in range(0, len(boundaryLayerPoints)):
    j = i*surfaceOrder
    boundaryLayerDividerLines.append(gmsh.model.geo.addLine(aerofoilPoints[j], boundaryLayerPoints[i]))
    gmsh.model.geo.mesh.setTransfiniteCurve(boundaryLayerDividerLines[i],boundaryLayerNumCells,"Progression",boundaryLayerProgression)

#Create boundary layer quadrilateral curve loops
boundaryLayerQuadLoops = []
for i in range(0, len(boundaryLayerPoints)-1):
    boundaryLayerQuadLoops.append(gmsh.model.geo.addCurveLoop([boundaryLayerLines[i],-boundaryLayerDividerLines[i+1],-aerofoilCurves[i],boundaryLayerDividerLines[i]]))

#Create boundary layer quadrilateral surfaces
boundaryLayerQuadSurfaces = []
for i in range(0, len(boundaryLayerQuadLoops)):
    boundaryLayerQuadSurfaces.append(gmsh.model.geo.addPlaneSurface([boundaryLayerQuadLoops[i]]))
    gmsh.model.geo.mesh.setTransfiniteSurface(boundaryLayerQuadSurfaces[i])
    gmsh.model.geo.mesh.setRecombine(2,boundaryLayerQuadSurfaces[i])

################################################################################
#   Wake Geometry
################################################################################

#Create points for wake
upperWakeStartPoint = boundaryLayerPoints[0]
lowerWakeStartPoint = boundaryLayerPoints[-1]
upperWakeEndPoint = gmsh.model.geo.addPoint(1+wakeLength, wakeEndThickness, 0, farFieldMeshSize)
lowerWakeEndPoint = gmsh.model.geo.addPoint(1+wakeLength, -wakeEndThickness, 0, farFieldMeshSize)

#Create lines for wake
upperWakeLine = gmsh.model.geo.addLine(upperWakeEndPoint, upperWakeStartPoint)
gmsh.model.geo.mesh.setTransfiniteCurve(upperWakeLine, wakeLengthPoints, "Progression", -1.05)
lowerWakeLine = gmsh.model.geo.addLine(lowerWakeStartPoint, lowerWakeEndPoint)
gmsh.model.geo.mesh.setTransfiniteCurve(lowerWakeLine, wakeLengthPoints, "Progression", 1.05)
rightWakeLine = gmsh.model.geo.addLine(lowerWakeEndPoint, upperWakeEndPoint)
gmsh.model.geo.mesh.setTransfiniteCurve(rightWakeLine, boundaryLayerNumCells*2)

#Create curveloop for wake
wakeLoop = gmsh.model.geo.addCurveLoop([upperWakeLine, -boundaryLayerDividerLines[0], trailingEdgeLine, boundaryLayerDividerLines[-1], lowerWakeLine, rightWakeLine])

#Create surface for wake
wakeSurface = gmsh.model.geo.addPlaneSurface([wakeLoop])
gmsh.model.geo.mesh.setTransfiniteSurface(wakeSurface, "Left", [upperWakeStartPoint, lowerWakeStartPoint, lowerWakeEndPoint, upperWakeEndPoint])
gmsh.model.geo.mesh.setRecombine(2, wakeSurface)

################################################################################
#   Far Field Geometry
################################################################################

#Create points for the far field
topLeft = gmsh.model.geo.addPoint(-farFieldSize/2+0.5, farFieldSize/2, 0, farFieldMeshSize)
topRight = gmsh.model.geo.addPoint(wakeLength+1, farFieldSize/2, 0, farFieldMeshSize)
bottomLeft = gmsh.model.geo.addPoint(-farFieldSize/2+0.5, -farFieldSize/2, 0, farFieldMeshSize)
bottomRight = gmsh.model.geo.addPoint(wakeLength+1, -farFieldSize/2, 0, farFieldMeshSize)

#Create lines for the far field
top = gmsh.model.geo.addLine(topRight, topLeft)
upperRight = gmsh.model.geo.addLine(upperWakeEndPoint, topRight)
lowerRight = gmsh.model.geo.addLine(bottomRight, lowerWakeEndPoint)
bottom = gmsh.model.geo.addLine(bottomLeft, bottomRight)
left = gmsh.model.geo.addLine(topLeft, bottomLeft)
farFieldLines = [top, left, bottom, lowerRight, -lowerWakeLine]
for i in range(0, len(boundaryLayerLines)):
    farFieldLines.append(-boundaryLayerLines[i])
farFieldLines.append(-upperWakeLine)
farFieldLines.append(upperRight)

#Create the curve loop for the far field
farFieldLoop = gmsh.model.geo.addCurveLoop(farFieldLines)

#Create the surface
farFieldSurface = gmsh.model.geo.addPlaneSurface([farFieldLoop])

#Synchronize CAD entities with the Gmsh model
gmsh.model.geo.synchronize()

################################################################################
#   Physical Groups
################################################################################

#Create physical groups
aerofoilGroup = gmsh.model.addPhysicalGroup(1, aerofoilCurves+[trailingEdgeLine])
gmsh.model.setPhysicalName(1, aerofoilGroup, "Aerofoil Boundary")
farFieldGroup = gmsh.model.addPhysicalGroup(1, [top, left, bottom, lowerRight, rightWakeLine, upperRight])
gmsh.model.setPhysicalName(1, farFieldGroup, "Far Field Boundary")
fluidGroup = gmsh.model.addPhysicalGroup(2, boundaryLayerQuadSurfaces + [wakeSurface, farFieldSurface])
gmsh.model.setPhysicalName(2, fluidGroup, "Fluid")

################################################################################
#   Output
################################################################################

#Set colours
# for i in range(0, len(boundaryLayerQuadSurfaces)):
#     gmsh.model.setColor([(2, boundaryLayerQuadSurfaces[i])], 255, 0, 0)
# gmsh.model.setColor([(2,volumeSurface)], 0, 255, 0)
# gmsh.option.setNumber("General.BackgroundGradient", 0)

#Generate and save mesh
gmsh.model.mesh.generate(2)
gmsh.write("NACA0012.msh")

#Visualize model
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

#Finalize
gmsh.finalize()