import naca_4_series_points
import gmsh
import sys
import aerofoil_scale

#Parameters
n = 25                          #Number of panels for each of the upper and lower surfaces of the aerofoil
farFieldSize = 10               #Width and height of the far field volume boundary
farFieldMeshSize = 1            #Target mesh size at the far field volume boundary
boundaryLayerThickness = 0.2    #Thickness of the boundary layer
boundaryLayerNumCells = 10      #Number of cells in the boundary layer in the direction normal to the aerofoil surface
boundaryLayerProgression = 1.5  #Growth rate of the cells in the boundary layer
wakeLength = 6
wakeStart = 1.5
wakeEndThickness = 2
wakeLengthPoints = 50
wakeHeightPoints = 20

################################################################################
#   Aerofoil Geometry
################################################################################

#Get the coordinates of the aerofoil
coords = naca_4_series_points.points(0, 0, 12, 100, "sin")

#Initialize Gmsh and name the model
gmsh.initialize()
gmsh.model.add("NACA0012")

#Create the aerofoil points
aerofoilPoints = []
trailingEdgePoint = gmsh.model.geo.addPoint(1, 0, 0)
for i in range(1, len(coords)-1):
    aerofoilPoints.append(gmsh.model.geo.addPoint(coords[i][0], coords[i][1], 0))

#Create the aerofoil lines
aerofoilSpline = gmsh.model.geo.addSpline(aerofoilPoints)
upperTrailingEdgeLine = gmsh.model.geo.addLine(trailingEdgePoint, aerofoilPoints[0])
lowerTrailingEdgeLine = gmsh.model.geo.addLine(aerofoilPoints[-1], trailingEdgePoint)

#Create the aerofoil curve loop
aerofoilLoop = gmsh.model.geo.addCurveLoop([upperTrailingEdgeLine, aerofoilSpline, lowerTrailingEdgeLine])

#Create boundary layer points
# leadingEdgeCircStartPoint = gmsh.model.geo.addPoint(0, 0.2, 0)
# leadingEdgeCircCentrePoint = gmsh.model.geo.addPoint(0, 0, 0)
# leadingEdgeCircEndPoint = gmsh.model.geo.addPoint(0, -0.2, 0)
# aboveTrailingEdgePoint = gmsh.model.geo.addPoint(1, 0.2, 0)
# belowTrailingEdgePoint = gmsh.model.geo.addPoint(1, -0.2, 0)

#Create the boundary layer lines
# upperBoundaryLayerLine = gmsh.model.geo.addLine(aboveTrailingEdgePoint, leadingEdgeCircStartPoint)
# gmsh.model.geo.mesh.setTransfiniteCurve(upperBoundaryLayerLine, 20)
# lowerBoundaryLayerLine = gmsh.model.geo.addLine(leadingEdgeCircEndPoint, belowTrailingEdgePoint)
# leadingEdgeCircLine = gmsh.model.geo.addCircleArc(leadingEdgeCircStartPoint, leadingEdgeCircCentrePoint, leadingEdgeCircEndPoint)
# upperRearBoundaryLayerLine = gmsh.model.geo.addLine(trailingEdgePoint, aboveTrailingEdgePoint)
# lowerRearBoundaryLayerLine = gmsh.model.geo.addLine(belowTrailingEdgePoint, trailingEdgePoint)

#Create the boundary layer curve loop
#boundaryLayerLoop = gmsh.model.geo.addCurveLoop([upperRearBoundaryLayerLine, upperBoundaryLayerLine, leadingEdgeCircLine, lowerBoundaryLayerLine, lowerRearBoundaryLayerLine])

#Create lines to divide up boundary layer into quadrilaterals
# boundaryLayerDividerLines = []
# for i in range(0, len(boundaryLayerPoints)):
#     boundaryLayerDividerLines.append(gmsh.model.geo.addLine(aerofoilPoints[i], boundaryLayerPoints[i]))
#     gmsh.model.geo.mesh.setTransfiniteCurve(boundaryLayerDividerLines[i],boundaryLayerNumCells,"Progression",boundaryLayerProgression)

#Create boundary layer quadrilateral curve loops
# boundaryLayerQuadLoops = []
# for i in range(0, len(boundaryLayerPoints)):
#     boundaryLayerQuadLoops.append(gmsh.model.geo.addCurveLoop([boundaryLayerLines[i],-boundaryLayerDividerLines[i],-aerofoilLines[i],boundaryLayerDividerLines[i-1]]))

#Create boundary layer quadrilateral surfaces
# boundaryLayerQuadSurfaces = []
# for i in range(0, len(boundaryLayerPoints)):
#     boundaryLayerQuadSurfaces.append(gmsh.model.geo.addPlaneSurface([boundaryLayerQuadLoops[i]]))
#     gmsh.model.geo.mesh.setTransfiniteSurface(boundaryLayerQuadSurfaces[i])
#     gmsh.model.geo.mesh.setRecombine(2,boundaryLayerQuadSurfaces[i])

################################################################################
#   Wake Geometry
################################################################################

#Create points for wake
upperWakeStartPoint = gmsh.model.geo.addPoint(wakeStart, boundaryLayerThickness, 0)
lowerWakeStartPoint = gmsh.model.geo.addPoint(wakeStart, -boundaryLayerThickness, 0)
upperWakeEndPoint = gmsh.model.geo.addPoint(wakeStart+wakeLength, wakeEndThickness, 0)
lowerWakeEndPoint = gmsh.model.geo.addPoint(wakeStart+wakeLength, -wakeEndThickness, 0)

#Create lines for wake
upperWakeLine = gmsh.model.geo.addLine(upperWakeEndPoint, upperWakeStartPoint)
gmsh.model.geo.mesh.setTransfiniteCurve(upperWakeLine, wakeLengthPoints)
lowerWakeLine = gmsh.model.geo.addLine(lowerWakeStartPoint, lowerWakeEndPoint)
gmsh.model.geo.mesh.setTransfiniteCurve(lowerWakeLine, wakeLengthPoints)
leftWakeLine = gmsh.model.geo.addLine(upperWakeStartPoint, lowerWakeStartPoint)
gmsh.model.geo.mesh.setTransfiniteCurve(leftWakeLine, wakeHeightPoints)
rightWakeLine = gmsh.model.geo.addLine(lowerWakeEndPoint, upperWakeEndPoint)
gmsh.model.geo.mesh.setTransfiniteCurve(rightWakeLine, wakeHeightPoints)

#Create curveloop for wake
wakeLoop = gmsh.model.geo.addCurveLoop([upperWakeLine, leftWakeLine, lowerWakeLine, rightWakeLine])

#Create surface for wake
wakeSurface = gmsh.model.geo.addPlaneSurface([wakeLoop])
gmsh.model.geo.mesh.setTransfiniteSurface(wakeSurface)
gmsh.model.geo.mesh.setRecombine(2,wakeSurface)

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