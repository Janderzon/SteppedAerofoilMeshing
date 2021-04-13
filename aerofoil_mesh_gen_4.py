import gmsh
import sys
import numpy as np
import math

################################################################################
#   Input Parameters
################################################################################

aerofoilName = "NACA0009"
blNumPoints = 30
Re = 100000
cMeshHeight = 1
alpha = 0

blThicknessProgression = 1.2
blTEProgression = 1.07
blLEBump = 6
oMeshProgression = 1.1
oCentreX = 0.5
wakeProgression = 1.2

################################################################################
#   Calculate Useful Parameters
################################################################################

#Calculate the the first cell height required to resolve laminar flow 
deltaLam = 1.227/(Re**0.75)
firstCellHeight = deltaLam/0.047
print("First cell height laminar flow: "+str(firstCellHeight))

#Calculate the the first cell height required to turbulent laminar flow 
deltaTurb = (13.1463**0.875)/(Re**0.9)
firstCellHeightTurb = deltaTurb/0.047
print("First cell height turbulent flow: "+str(firstCellHeightTurb))

#Calculate the number of points required in the normal direction in the C mesh to acheive the desired first cell height
cNumPoints = round(math.log(1-(cMeshHeight*(1-oMeshProgression)/firstCellHeight))/math.log(oMeshProgression))
print(math.log(1-(cMeshHeight*(1-oMeshProgression)/firstCellHeight))/math.log(oMeshProgression))
cMeshHeight = firstCellHeight*(1-oMeshProgression**cNumPoints)/(1-oMeshProgression)
print(cMeshHeight)

#Calculate the far field cell size
ffCellSize = firstCellHeight*oMeshProgression**(cNumPoints-1)

#Calculate the number of points required in the wake for a smooth transition from trailing edge to wake
teCellWidth = 0.7*(1-blTEProgression)/(1-blTEProgression**(blNumPoints*0.7-1))
wakeNumPoints = round(math.log(ffCellSize/teCellWidth)/math.log(wakeProgression)+1)
wakeProgressionLength = teCellWidth*(1-wakeProgression**wakeNumPoints)/(1-wakeProgression)
ffNumPoints = round((2*cMeshHeight-wakeProgressionLength)/ffCellSize)


#Density of far field points
ffPointsPerLength = 1/(firstCellHeight*cMeshHeight*40)
ffPointsPerLength = 5/cMeshHeight

#Number of points for the leading and trailing edge sections over the aerofoil
blTENumPoints = round((1-0.35)*blNumPoints)
blLENumPoints = round(0.35*blNumPoints)

################################################################################
#   Initialize Gmsh
################################################################################

#Initialize Gmsh and name the model
gmsh.initialize()
gmsh.model.add(aerofoilName)

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
#   Reference Geometry
################################################################################

origin = gmsh.model.geo.addPoint(0, 0, 0)
trailingEdge = gmsh.model.geo.addPoint(1, 0, 0)
thickestX = gmsh.model.geo.addPoint(thickestXCoord, 0, 0)
oCentre = gmsh.model.geo.addPoint(oCentreX, 0, 0)

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

gmsh.model.geo.rotate([(1, teUpperSpline), (1, leSpline), (1, teLowerSpline)], thickestXCoord, 0, 0, 0, 0, 1, math.radians(-alpha))

#Set aerofoil transfinite curves
gmsh.model.geo.mesh.setTransfiniteCurve(teUpperSpline, blTENumPoints, "Progression", blTEProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(leSpline, blLENumPoints*2+1, "Bump", blLEBump)
gmsh.model.geo.mesh.setTransfiniteCurve(teLowerSpline, blTENumPoints, "Progression", -blTEProgression)

################################################################################
#   C Mesh Geometry
################################################################################

#Create C mesh points
oTopPoint = gmsh.model.geo.addPoint(thickestXCoord, cMeshHeight, 0)
oLeftPoint = gmsh.model.geo.addPoint(-cMeshHeight+thickestXCoord, 0, 0)
oBottomPoint = gmsh.model.geo.addPoint(thickestXCoord, -cMeshHeight, 0)
oRightPoint = gmsh.model.geo.addPoint(cMeshHeight+1, 0, 0)

#Create C mesh lines
oTopLeftArc = gmsh.model.geo.addCircleArc(oTopPoint, thickestX, oLeftPoint)
oBottomLeftArc = gmsh.model.geo.addCircleArc(oLeftPoint, thickestX, oBottomPoint)
oBottomRightArc = gmsh.model.geo.addEllipseArc(oBottomPoint, thickestX, oRightPoint, oRightPoint)
oTopRightArc = gmsh.model.geo.addEllipseArc(oRightPoint, thickestX, oRightPoint, oTopPoint)
oTopDivideLine = gmsh.model.geo.addLine(lePoints[0], oTopPoint)
oBottomDivideLine = gmsh.model.geo.addLine(oBottomPoint, teLowerPoints[0])

#Set transfinite C mesh lines
gmsh.model.geo.mesh.setTransfiniteCurve(oTopLeftArc, blLENumPoints+1)
gmsh.model.geo.mesh.setTransfiniteCurve(oBottomLeftArc, blLENumPoints+1)
gmsh.model.geo.mesh.setTransfiniteCurve(oBottomRightArc, blTENumPoints)
gmsh.model.geo.mesh.setTransfiniteCurve(oTopRightArc, blTENumPoints)
gmsh.model.geo.mesh.setTransfiniteCurve(oTopDivideLine, cNumPoints, "Progression", oMeshProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(oBottomDivideLine, cNumPoints, "Progression", -oMeshProgression)
# gmsh.model.geo.mesh.setTransfiniteCurve(cMidLeftLine, wakeNumPoints, "Progression", wakeProgression)
# gmsh.model.geo.mesh.setTransfiniteCurve(cMidRightLine, ffNumPoints)
# gmsh.model.geo.mesh.setTransfiniteCurve(cTopLine, wakeNumPoints+ffNumPoints-1)
# gmsh.model.geo.mesh.setTransfiniteCurve(cBottomLine, wakeNumPoints+ffNumPoints-1)
# gmsh.model.geo.mesh.setTransfiniteCurve(cTopLeftLine, blTENumPoints, "Progression", blTEProgression)
# gmsh.model.geo.mesh.setTransfiniteCurve(cBottomLeftLine, blTENumPoints, "Progression", -blTEProgression)
# gmsh.model.geo.mesh.setTransfiniteCurve(cTopRightLine, cNumPoints, "Progression", oMeshProgression)
# gmsh.model.geo.mesh.setTransfiniteCurve(cBottomRightLine, cNumPoints, "Progression", -oMeshProgression)

#Create C mesh curve loop
oLeftLoop = gmsh.model.geo.addCurveLoop([oTopLeftArc, oBottomLeftArc, oBottomDivideLine, -leSpline, oTopDivideLine])
oRightLoop = gmsh.model.geo.addCurveLoop([oTopRightArc, -oTopDivideLine, -teUpperSpline, -teLowerSpline, -oBottomDivideLine, oBottomRightArc])

#Create C mesh surface
oLeftSurface = gmsh.model.geo.addPlaneSurface([oLeftLoop])
oRightSurface = gmsh.model.geo.addPlaneSurface([oRightLoop])
gmsh.model.geo.mesh.setTransfiniteSurface(oLeftSurface, "Left", [oTopPoint, lePoints[0], teLowerPoints[0], oBottomPoint])
gmsh.model.geo.mesh.setTransfiniteSurface(oRightSurface, "Left", [oTopPoint, lePoints[0], teLowerPoints[0], oBottomPoint])
gmsh.model.geo.mesh.setRecombine(2, oLeftSurface)
gmsh.model.geo.mesh.setRecombine(2, oRightSurface)

################################################################################
#   Far Field Geometry
################################################################################

# #Create points for the far field
# ffTopRightPoint = gmsh.model.geo.addPoint(2*ffBoundaryScaling, ffBoundaryScaling, 0)
# ffTopLeftPoint = gmsh.model.geo.addPoint(-ffBoundaryScaling, ffBoundaryScaling, 0)
# ffBottomLeftPoint = gmsh.model.geo.addPoint(-ffBoundaryScaling, -ffBoundaryScaling, 0)
# ffBottomRightPoint = gmsh.model.geo.addPoint(2*ffBoundaryScaling, -ffBoundaryScaling, 0)

# #Create lines for the far field
# ffLeftLine = gmsh.model.geo.addLine(ffTopLeftPoint, ffBottomLeftPoint)
# ffBottomLine = gmsh.model.geo.addLine(ffBottomLeftPoint, ffBottomRightPoint)
# ffTopLine = gmsh.model.geo.addLine(ffTopRightPoint, ffTopLeftPoint)
# ffTopRightLine = gmsh.model.geo.addLine(cTopRightPoint, ffTopRightPoint)
# ffBottomRightLine = gmsh.model.geo.addLine(ffBottomRightPoint, cBottomRightPoint)

# #Set transfinite curves for far field
# gmsh.model.geo.mesh.setTransfiniteCurve(ffLeftLine, round(ffPointsPerLength*ffBoundaryScaling*2))
# gmsh.model.geo.mesh.setTransfiniteCurve(ffBottomLine, round(ffPointsPerLength*ffBoundaryScaling*3))
# gmsh.model.geo.mesh.setTransfiniteCurve(ffTopLine, round(ffPointsPerLength*ffBoundaryScaling*3))
# gmsh.model.geo.mesh.setTransfiniteCurve(ffTopRightLine, round(ffPointsPerLength*ffBoundaryScaling))
# gmsh.model.geo.mesh.setTransfiniteCurve(ffBottomRightLine, round(ffPointsPerLength*ffBoundaryScaling))

# #Create the curve loop for the far field
# farFieldLoop = gmsh.model.geo.addCurveLoop([ffTopLine, ffLeftLine, ffBottomLine, ffBottomRightLine, -cBottomLine, -cBottomLeftLine, -cBottomArc, -cTopArc, -cTopLeftLine, -cTopLine, ffTopRightLine])

# #Create the surface
# farFieldSurface = gmsh.model.geo.addPlaneSurface([farFieldLoop])

################################################################################
#   Synchronize
################################################################################

#Synchronize CAD entities with the Gmsh model
gmsh.model.geo.synchronize()

################################################################################
#   Physical Groups
################################################################################

# #Create physical groups
# aerofoilGroup = gmsh.model.addPhysicalGroup(1, [teUpperSpline, leSpline, teLowerSpline])
# gmsh.model.setPhysicalName(1, aerofoilGroup, "aerofoilwall")
# farFieldGroup = gmsh.model.addPhysicalGroup(1, [wakeRightLine, ffTopLine, ffLeftLine, ffBottomLine, ffBottomRightLine, ffTopLine2, ffRightLine2])
# gmsh.model.setPhysicalName(1, farFieldGroup, "farfieldboundary")
# fluidGroup = gmsh.model.addPhysicalGroup(2, [blLESurface, blTopTESurface, blBottomTESurface, wakeSurface, srSurface, farFieldSurface, farFieldSurface2])
# gmsh.model.setPhysicalName(2, fluidGroup, "fluid")

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
gmsh.option.setNumber("Mesh.HighOrderOptimize", 4)
gmsh.option.setNumber("Mesh.NumSubEdges", 10)
gmsh.model.mesh.generate(2)
gmsh.write(aerofoilName+".msh")

#Visualize model
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

#Finalize
gmsh.finalize()