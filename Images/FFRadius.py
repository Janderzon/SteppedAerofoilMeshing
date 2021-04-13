import gmsh
import sys
import numpy as np
import math

################################################################################
#   Input Parameters
################################################################################

aerofoilName = "NACA0009"
Re = 80000
alpha = 0

ffBoundaryScaling = 2
blNumPoints = 100
cNoseNearWallHeight = 0.001
outerLayerHeight = 1.25
outerLayerLength = 11
outerLayerCellSizeScaling = 6

blTEProgression = 1.04
blLEBump = 10
cProgression = 1.2
wakeProgression = 1.2
blPointsSplit = 0.45
cMeshHeight = 0.175
cWakeLength = 1

ffDensity = 20

################################################################################
#   Aerofoil Coordinates
################################################################################

#Get the coordinates of the aerofoil
aerofoilData = np.loadtxt(aerofoilName+".dat")
coords = aerofoilData

#Find the points of maximum thickness on the aerofoil.
thickestIndex = 0
for i in range(0, len(coords)):
    if(coords[i][1]>coords[thickestIndex][1]):
        thickestIndex = i
thickestXCoord = coords[thickestIndex][0]
thickestYCoord = coords[thickestIndex][1]

################################################################################
#   Calculate Useful Parameters
################################################################################

#Number of points for the leading and trailing edge sections over the aerofoil
blTENumPoints = round((1-blPointsSplit)*blNumPoints)
blLENumPoints = round(blPointsSplit*blNumPoints)

#Calculate the the first cell height required to resolve laminar flow 
deltaLam = 1.3016/(Re**0.75)
firstCellHeight = deltaLam/0.047
print("First cell height laminar flow: "+str(firstCellHeight))

cNoseNearWallHeight = firstCellHeight*0.5

#Calculate the the first cell height required to resolve turbulent flow 
deltaTurb = (13.1463**0.875)/(Re**0.9)
firstCellHeightTurb = deltaTurb/0.047
print("First cell height turbulent flow: "+str(firstCellHeightTurb))

#Calculate the number of points required in the normal direction in the C mesh to acheive the desired first cell height
cNumPoints = round(math.log(1-(cMeshHeight*(1-cProgression)/firstCellHeight))/math.log(cProgression))+1
cMeshHeight = firstCellHeight*(1-cProgression**(cNumPoints-1))/(1-cProgression)

#Calculate the far field cell size
cFarCellSize = firstCellHeight*cProgression**(cNumPoints-2)

#Calculate the number of points required in the wake for a smooth transition from trailing edge to wake
teCellWidth = (1-thickestXCoord)*(1-blTEProgression)/(1-blTEProgression**(blTENumPoints-1))
wakeNumPoints = round(math.log(cFarCellSize/teCellWidth)/math.log(wakeProgression))+2
wakeProgressionLength = teCellWidth*(1-wakeProgression**(wakeNumPoints-1))/(1-wakeProgression)

#Distance of C mesh nose from leading edge
cNoseDistance = cNoseNearWallHeight*(1-cProgression**(cNumPoints-1))/(1-cProgression)

################################################################################
#   Initialize Gmsh
################################################################################

#Initialize Gmsh and name the model
gmsh.initialize()
gmsh.model.add(aerofoilName)

################################################################################
#   Reference Geometry
################################################################################

origin = gmsh.model.geo.addPoint(0, 0, 0)
trailingEdge = gmsh.model.geo.addPoint(1, 0, 0)
thickestX = gmsh.model.geo.addPoint(thickestXCoord, 0, 0)

################################################################################
#   Aerofoil Geometry
################################################################################

#Create the aerofoil points
aerofoilPoints = []
for i in range(0, len(coords)):
    aerofoilPoints.append(gmsh.model.geo.addPoint(coords[i][0], coords[i][1], 0))

#Divide up the aerofoil point into leading and upper and lower trailing edges
lePoints = []
teUpperPoints = []
teLowerPoints = []
teUpperPoints.append(trailingEdge)
for i in range(1, thickestIndex+1):
    teUpperPoints.append(aerofoilPoints[i])
for i in range(thickestIndex, len(coords)-thickestIndex):
    lePoints.append(aerofoilPoints[i])
for i in range(len(coords)-thickestIndex-1, len(coords)-1):
    teLowerPoints.append(aerofoilPoints[i])
teLowerPoints.append(trailingEdge)

#Create the aerofoil lines
teUpperSpline = gmsh.model.geo.addSpline(teUpperPoints)
leSpline = gmsh.model.geo.addSpline(lePoints)
teLowerSpline = gmsh.model.geo.addSpline(teLowerPoints)

#Set aerofoil transfinite curves
gmsh.model.geo.mesh.setTransfiniteCurve(teUpperSpline, blTENumPoints, "Progression", blTEProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(leSpline, blLENumPoints*2+1, "Bump", blLEBump)
gmsh.model.geo.mesh.setTransfiniteCurve(teLowerSpline, blTENumPoints, "Progression", -blTEProgression)

################################################################################
#   Far Field Geometry
################################################################################

#Create points for the far field
ffTopPoint = gmsh.model.geo.addPoint(1, ffBoundaryScaling, 0, ffDensity)
ffLeftPoint = gmsh.model.geo.addPoint(1-ffBoundaryScaling, 0, 0, ffDensity)
ffBottomPoint = gmsh.model.geo.addPoint(1, -ffBoundaryScaling, 0, ffDensity)
ffRightPoint = gmsh.model.geo.addPoint(1+ffBoundaryScaling, 0, 0, ffDensity)
edgePoint = gmsh.model.geo.addPoint(1+ffBoundaryScaling/math.sqrt(2), ffBoundaryScaling/math.sqrt(2), 0)
edgePoint1 = gmsh.model.geo.addPoint(1+ffBoundaryScaling/math.sqrt(2), ffBoundaryScaling/math.sqrt(2)-0.05, 0)
edgePoint2 = gmsh.model.geo.addPoint(0.95+ffBoundaryScaling/math.sqrt(2), ffBoundaryScaling/math.sqrt(2), 0)

#Create lines for the far field
ffTopLeftArc = gmsh.model.geo.addCircleArc(ffTopPoint, trailingEdge, ffLeftPoint)
ffBottomLeftArc = gmsh.model.geo.addCircleArc(ffLeftPoint, trailingEdge, ffBottomPoint)
ffTopRightArc = gmsh.model.geo.addCircleArc(ffBottomPoint, trailingEdge, ffRightPoint)
ffBottomRightArc = gmsh.model.geo.addCircleArc(ffRightPoint, trailingEdge, ffTopPoint)
gmsh.model.geo.addLine(aerofoilPoints[0], edgePoint)
gmsh.model.geo.addLine(edgePoint1, edgePoint)
gmsh.model.geo.addLine(edgePoint2, edgePoint)

#Create the curve loop for the far field
farFieldLoop = gmsh.model.geo.addCurveLoop([ffTopLeftArc, ffBottomLeftArc, ffBottomRightArc, ffTopRightArc])

################################################################################
#   Synchronize
################################################################################

#Synchronize CAD entities with the Gmsh model
gmsh.model.geo.synchronize()

################################################################################
#   Output
################################################################################

#Set colours
# for i in range(0, len(boundaryLayerQuadSurfaces)):
#     gmsh.model.setColor([(2, boundaryLayerQuadSurfaces[i])], 255, 0, 0)
# gmsh.model.setColor([(2,volumeSurface)], 0, 255, 0)
# gmsh.option.setNumber("General.BackgroundGradient", 0)

#Generate and save mesh
gmsh.option.setNumber("Mesh.ElementOrder", 4)
#gmsh.option.setNumber("Mesh.HighOrderOptimize", 4)
gmsh.option.setNumber("Mesh.NumSubEdges", 10)
gmsh.model.mesh.generate(2)
gmsh.write(aerofoilName+".msh")

#Visualize model
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

#Finalize
gmsh.finalize()