import naca_4_series_points
import gmsh
import sys
import math
import aerofoil_scale

#Parameters
n = 50                          #Number points on the aerofoil surface
surfaceOrder = 2
farFieldSize = 10               #Width and height of the far field volume boundary
farFieldMeshSize = 1            #Target mesh size at the far field volume boundary
boundaryLayerThickness = 0.2    #Thickness of the boundary layer
boundaryLayerNumCells = 10      #Number of cells in the boundary layer in the direction normal to the aerofoil surface
boundaryLayerProgression = 1.3  #Growth rate of the cells in the boundary layer
wakeLength = 6
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
    print(angle)
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
upperWakeEndPoint = gmsh.model.geo.addPoint(1+wakeLength, wakeEndThickness, 0)
lowerWakeEndPoint = gmsh.model.geo.addPoint(1+wakeLength, -wakeEndThickness, 0)

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

#Create points for the volume boundary
# topLeft = gmsh.model.geo.addPoint(-farFieldSize/2+0.5, farFieldSize/2, 0, farFieldMeshSize)
# topRight = gmsh.model.geo.addPoint(farFieldSize/2+0.5, farFieldSize/2, 0, farFieldMeshSize)
# bottomLeft = gmsh.model.geo.addPoint(-farFieldSize/2+0.5, -farFieldSize/2, 0, farFieldMeshSize)
# bottomRight = gmsh.model.geo.addPoint(farFieldSize/2+0.5, -farFieldSize/2, 0, farFieldMeshSize)

#Create lines for the volume boundary
# top = gmsh.model.geo.addLine(topLeft, topRight)
# right = gmsh.model.geo.addLine(topRight, bottomRight)
# bottom = gmsh.model.geo.addLine(bottomRight, bottomLeft)
# left = gmsh.model.geo.addLine(bottomLeft, topLeft)

#Create the curve loop for the volume boundary
# volumeLoop = gmsh.model.geo.addCurveLoop([top, right, bottom, left])

#Create the surface
#volumeSurface = gmsh.model.geo.addPlaneSurface([boundaryLayerLoop, volumeLoop])
#volumeSurface = gmsh.model.geo.addPlaneSurface([volumeLoop,boundaryLayerLoop])

#Synchronize CAD entities with the Gmsh model
gmsh.model.geo.synchronize()

#Create physical groups
# aerofoilGroup = gmsh.model.addPhysicalGroup(1, [aerofoilLoop])
# gmsh.model.setPhysicalName(1, aerofoilGroup, "Aerofoil Surface")
# volumeGroup = gmsh.model.addPhysicalGroup(1, [volumeLoop])
# gmsh.model.setPhysicalName(1, volumeGroup, "Volume Boundary")

#Set colours
# for i in range(0, len(boundaryLayerQuadSurfaces)):
#     gmsh.model.setColor([(2, boundaryLayerQuadSurfaces[i])], 255, 0, 0)
# gmsh.model.setColor([(2,volumeSurface)], 0, 255, 0)
# gmsh.option.setNumber("General.BackgroundGradient", 0)

#Generate and save mesh
gmsh.model.mesh.generate(2)
#gmsh.write("NACA0012.msh")

#Visualize model
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

#Finalize
gmsh.finalize()