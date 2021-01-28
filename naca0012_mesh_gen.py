import naca_4_series_points
import gmsh
import sys
import math

################################################################################
#   Input Parameters
################################################################################

n = 100                         #Number points on the aerofoil surface
surfaceOrder = 1
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
coords = naca_4_series_points.points(0, 0, 12, n*surfaceOrder, "sin")

#Initialize Gmsh and name the model
gmsh.initialize()
gmsh.model.add("NACA0012")

#Create the aerofoil points
aerofoilPoints = []
for i in range(0, len(coords)):
    aerofoilPoints.append(gmsh.model.geo.addPoint(coords[i][0], coords[i][1], 0))

#Divide up the aerofoil point into leading and upper and lower trailing edges
thickestIndex = 0
for i in range(0, len(coords)):
    if(coords[i][1]>coords[thickestIndex][1]):
        thickestIndex = i
lePoints = []
teUpperPoints = []
teLowerPoints = []
for i in range(0, thickestIndex+1):
    teUpperPoints.append(aerofoilPoints[i])
for i in range(thickestIndex, len(coords)-thickestIndex):
    lePoints.append(aerofoilPoints[i])
for i in range(len(coords)-thickestIndex-1, len(coords)):
    teLowerPoints.append(aerofoilPoints[i])

#Create the aerofoil lines
teLine = gmsh.model.geo.addLine(aerofoilPoints[-1], aerofoilPoints[0])
teUpperSpline = gmsh.model.geo.addSpline(teUpperPoints)
leSpline = gmsh.model.geo.addSpline(lePoints)
teLowerSpline = gmsh.model.geo.addSpline(teLowerPoints)
gmsh.model.geo.mesh.setTransfiniteCurve(teUpperSpline, 100, "Progression", 1.01)
gmsh.model.geo.mesh.setTransfiniteCurve(leSpline, 100, "Bump", 1.5)
gmsh.model.geo.mesh.setTransfiniteCurve(teLowerSpline, 100, "Progression", -1.01)

#Create the aerofoil curve loop
aerofoilLoop = gmsh.model.geo.addCurveLoop([teUpperSpline, leSpline, teLowerSpline, teLine])

################################################################################
#   Boundary Layer Geometry
################################################################################

#Create boundary layer points
blTopRightPoint = gmsh.model.geo.addPoint(1, boundaryLayerThickness, 0)
blBottomRightPoint = gmsh.model.geo.addPoint(1, -boundaryLayerThickness, 0)
blTopMidPoint = gmsh.model.geo.addPoint(coords[thickestIndex][0], boundaryLayerThickness, 0)
blBottomMidPoint = gmsh.model.geo.addPoint(coords[thickestIndex][0], -boundaryLayerThickness, 0)

#Create the boundary layer lines
blTopTELine = gmsh.model.geo.addLine(blTopRightPoint, blTopMidPoint)
blBottomTELine = gmsh.model.geo.addLine(blBottomMidPoint, blBottomRightPoint)
blTopRightLine = gmsh.model.geo.addLine(aerofoilPoints[0], blTopRightPoint)
blBottomRightLine = gmsh.model.geo.addLine(blBottomRightPoint, aerofoilPoints[-1])
blTopMidLine = gmsh.model.geo.addLine(blTopMidPoint, aerofoilPoints[thickestIndex])
blBottomMidLine = gmsh.model.geo.addLine(aerofoilPoints[len(coords)-thickestIndex-1], blBottomMidPoint)
gmsh.model.geo.mesh.setTransfiniteCurve(blTopTELine, 100)
gmsh.model.geo.mesh.setTransfiniteCurve(blBottomTELine, 100)
gmsh.model.geo.mesh.setTransfiniteCurve(blTopRightLine, boundaryLayerNumCells, "Progression", boundaryLayerProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(blBottomRightLine, boundaryLayerNumCells, "Progression", -boundaryLayerProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(blTopMidLine, boundaryLayerNumCells, "Progression", -boundaryLayerProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(blBottomMidLine, boundaryLayerNumCells, "Progression", boundaryLayerProgression)

#Create boundary layer curve loop
blTopTELoop = gmsh.model.geo.addCurveLoop([blTopTELine, blTopMidLine, -teUpperSpline, blTopRightLine])
blBottomTELoop = gmsh.model.geo.addCurveLoop([-teLowerSpline, blBottomMidLine, blBottomTELine, blBottomRightLine])

#Create boundary layer surfaces
blTopTESurface = gmsh.model.geo.addPlaneSurface([blTopTELoop])
blBottomTESurface = gmsh.model.geo.addPlaneSurface([blBottomTELoop])
gmsh.model.geo.mesh.setTransfiniteSurface(blTopTESurface)
gmsh.model.geo.mesh.setTransfiniteSurface(blBottomTESurface)
gmsh.model.geo.mesh.setRecombine(2, blTopTESurface)
gmsh.model.geo.mesh.setRecombine(2, blBottomTESurface)

################################################################################
#   Wake Geometry
################################################################################

# #Create points for wake
# upperWakeStartPoint = boundaryLayerPoints[0]
# lowerWakeStartPoint = boundaryLayerPoints[-1]
# upperWakeEndPoint = gmsh.model.geo.addPoint(1+wakeLength, wakeEndThickness, 0, farFieldMeshSize)
# lowerWakeEndPoint = gmsh.model.geo.addPoint(1+wakeLength, -wakeEndThickness, 0, farFieldMeshSize)

# #Create lines for wake
# upperWakeLine = gmsh.model.geo.addLine(upperWakeEndPoint, upperWakeStartPoint)
# gmsh.model.geo.mesh.setTransfiniteCurve(upperWakeLine, wakeLengthPoints, "Progression", -1.05)
# lowerWakeLine = gmsh.model.geo.addLine(lowerWakeStartPoint, lowerWakeEndPoint)
# gmsh.model.geo.mesh.setTransfiniteCurve(lowerWakeLine, wakeLengthPoints, "Progression", 1.05)
# rightWakeLine = gmsh.model.geo.addLine(lowerWakeEndPoint, upperWakeEndPoint)
# gmsh.model.geo.mesh.setTransfiniteCurve(rightWakeLine, boundaryLayerNumCells*2)

# #Create curveloop for wake
# wakeLoop = gmsh.model.geo.addCurveLoop([upperWakeLine, -boundaryLayerDividerLines[0], trailingEdgeLine, boundaryLayerDividerLines[-1], lowerWakeLine, rightWakeLine])

# #Create surface for wake
# wakeSurface = gmsh.model.geo.addPlaneSurface([wakeLoop])
# gmsh.model.geo.mesh.setTransfiniteSurface(wakeSurface, "Left", [upperWakeStartPoint, lowerWakeStartPoint, lowerWakeEndPoint, upperWakeEndPoint])
# gmsh.model.geo.mesh.setRecombine(2, wakeSurface)

################################################################################
#   Far Field Geometry
################################################################################

#Create points for the far field
# topLeft = gmsh.model.geo.addPoint(-farFieldSize/2+0.5, farFieldSize/2, 0, farFieldMeshSize)
# topRight = gmsh.model.geo.addPoint(wakeLength+1, farFieldSize/2, 0, farFieldMeshSize)
# bottomLeft = gmsh.model.geo.addPoint(-farFieldSize/2+0.5, -farFieldSize/2, 0, farFieldMeshSize)
# bottomRight = gmsh.model.geo.addPoint(wakeLength+1, -farFieldSize/2, 0, farFieldMeshSize)

#Create lines for the far field
#op = gmsh.model.geo.addLine(topRight, topLeft)
#upperRight = gmsh.model.geo.addLine(upperWakeEndPoint, topRight)
#lowerRight = gmsh.model.geo.addLine(bottomRight, lowerWakeEndPoint)
# right = gmsh.model.geo.addLine(bottomRight, topRight)
# bottom = gmsh.model.geo.addLine(bottomLeft, bottomRight)
# left = gmsh.model.geo.addLine(topLeft, bottomLeft)
# farFieldLines = [top, left, bottom, right]
# for i in range(0, len(boundaryLayerLines)):
#     farFieldLines.append(-boundaryLayerLines[i])
# farFieldLines.append(-upperWakeLine)
# farFieldLines.append(upperRight)

#Create the curve loop for the far field
#farFieldLoop = gmsh.model.geo.addCurveLoop(farFieldLines)

#Create the surface
#farFieldSurface = gmsh.model.geo.addPlaneSurface([farFieldLoop, aerofoilLoop])

#Synchronize CAD entities with the Gmsh model
gmsh.model.geo.synchronize()

################################################################################
#   Physical Groups
################################################################################

#Create physical groups
#aerofoilGroup = gmsh.model.addPhysicalGroup(1, [upperAerofoilSpline, lowerAerofoilSpline, trailingEdgeLine])
#gmsh.model.setPhysicalName(1, aerofoilGroup, "aerofoilboundary")
#farFieldGroup = gmsh.model.addPhysicalGroup(1, farFieldLines)
#gmsh.model.setPhysicalName(1, farFieldGroup, "farfieldboundary")
#fluidGroup = gmsh.model.addPhysicalGroup(2, [farFieldSurface])
#gmsh.model.setPhysicalName(2, fluidGroup, "fluid")

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