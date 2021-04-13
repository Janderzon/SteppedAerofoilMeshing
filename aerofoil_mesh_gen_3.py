import gmsh
import sys
import numpy as np
import math

################################################################################
#   Input Parameters
################################################################################

aerofoilName = "NACA0009"
blNumPoints = 30
ffBoundaryScaling = 10
Re = 100000
cMeshHeight = 10
alpha = 0

blThicknessProgression = 1.2
blTEProgression = 1.07
blLEBump = 6
ffProgression = 1.2
cMeshProgression = 1.1
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
cNumPoints = round(math.log(1-(cMeshHeight*(1-cMeshProgression)/firstCellHeight))/math.log(cMeshProgression))

#Calculate the number of points required in the wake for a smooth transition from trailing edge to wake
teCellWidth = 0.7*(1-blTEProgression)/(1-blTEProgression**(blNumPoints*0.7-1))
wakeNumPoints = round(math.log(1-(2*ffBoundaryScaling*(1-wakeProgression)/teCellWidth))/math.log(wakeProgression))

#Density of far field points
ffPointsPerLength = 1/(firstCellHeight*ffBoundaryScaling*40)
ffPointsPerLength = 5/ffBoundaryScaling

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
#gmsh.model.geo.translate([(1, teUpperSpline), (1, leSpline), (1, teLowerSpline)], 0, -0.5, 0)

#Set aerofoil transfinite curves
gmsh.model.geo.mesh.setTransfiniteCurve(teUpperSpline, blTENumPoints, "Progression", blTEProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(leSpline, blLENumPoints*2+1, "Bump", blLEBump)
gmsh.model.geo.mesh.setTransfiniteCurve(teLowerSpline, blTENumPoints, "Progression", -blTEProgression)

################################################################################
#   C Mesh Geometry
################################################################################

#Create C mesh points
cTopLeftPoint = gmsh.model.geo.addPoint(thickestXCoord, cMeshHeight, 0)
cLeftPoint = gmsh.model.geo.addPoint(-cMeshHeight*math.cos(math.radians(alpha))+thickestXCoord, cMeshHeight*math.sin(math.radians(alpha)), 0)
cBottomLeftPoint = gmsh.model.geo.addPoint(-cMeshHeight*math.sin(math.radians(alpha))+thickestXCoord, -cMeshHeight*math.cos(math.radians(alpha)), 0)
cTopRightPoint = gmsh.model.geo.addPoint(2*ffBoundaryScaling, cMeshHeight, 0)
cRightPoint = gmsh.model.geo.addPoint(2*ffBoundaryScaling, -0.5*ffBoundaryScaling*math.sin(math.radians(alpha)), 0)
cBottomRightPoint = gmsh.model.geo.addPoint(2*ffBoundaryScaling, -cMeshHeight-ffBoundaryScaling*math.sin(math.radians(alpha)), 0)
cMidPoint = gmsh.model.geo.addPoint(ffBoundaryScaling, 0, 0)

#Create C mesh lines
cTopArc = gmsh.model.geo.addCircleArc(cTopLeftPoint, thickestX, cLeftPoint)
cBottomArc = gmsh.model.geo.addCircleArc(cLeftPoint, thickestX, cBottomLeftPoint)
cTopMidLine = gmsh.model.geo.addLine(lePoints[0], cTopLeftPoint)
cBottomMidLine = gmsh.model.geo.addLine(cBottomLeftPoint, teLowerPoints[0])
cTopLine = gmsh.model.geo.addLine(cTopRightPoint, cTopLeftPoint)
cBottomLine = gmsh.model.geo.addLine(cBottomLeftPoint, cBottomRightPoint)
cTopRightLine = gmsh.model.geo.addLine(cRightPoint, cTopRightPoint)
cBottomRightLine = gmsh.model.geo.addLine(cBottomRightPoint, cRightPoint)
cMidRightLine = gmsh.model.geo.addLine(cMidPoint, cRightPoint)
cMidLeftLine = gmsh.model.geo.addLine(trailingEdge, cMidPoint)

#Set transfinite C mesh lines
gmsh.model.geo.mesh.setTransfiniteCurve(cTopArc, blLENumPoints+1)
gmsh.model.geo.mesh.setTransfiniteCurve(cBottomArc, blLENumPoints+1)
gmsh.model.geo.mesh.setTransfiniteCurve(cMidLeftLine, round(wakeNumPoints/2), "Progression", wakeProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(cMidRightLine, wakeNumPoints-round(wakeNumPoints/2))
gmsh.model.geo.mesh.setTransfiniteCurve(cTopLine, wakeNumPoints+blTENumPoints-2)
gmsh.model.geo.mesh.setTransfiniteCurve(cBottomLine, wakeNumPoints+blTENumPoints)
gmsh.model.geo.mesh.setTransfiniteCurve(cTopRightLine, cNumPoints)
gmsh.model.geo.mesh.setTransfiniteCurve(cBottomRightLine, cNumPoints)
gmsh.model.geo.mesh.setTransfiniteCurve(cTopMidLine, cNumPoints, "Progression", cMeshProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(cBottomMidLine, cNumPoints, "Progression", -cMeshProgression)

#Create C mesh curve loop
# cLoop = gmsh.model.geo.addCurveLoop([cTopLine, cTopLeftLine, cTopArc, cBottomArc, cBottomLeftLine, cBottomLine, cBottomRightLine, -cMidLine, -teLowerSpline, -leSpline, -teUpperSpline, cMidLine, cTopRightLine])
cNoseLoop = gmsh.model.geo.addCurveLoop([cTopArc, cBottomArc, cBottomMidLine, -leSpline, cTopMidLine])
cTopLoop = gmsh.model.geo.addCurveLoop([cTopLine, -cTopMidLine, -teUpperSpline, cMidLeftLine, cMidRightLine, cTopRightLine])
cBottomLoop = gmsh.model.geo.addCurveLoop([-cMidRightLine, -cMidLeftLine, -teLowerSpline, -cBottomMidLine, cBottomLine, cBottomRightLine])


# #Create C mesh surface
# cSurface = gmsh.model.geo.addPlaneSurface([cLoop])
# gmsh.model.geo.mesh.setTransfiniteSurface(cSurface, "Left", [cRightPoint, cTopRightPoint, cBottomRightPoint, cRightPoint])
# gmsh.model.geo.mesh.setRecombine(2, cSurface)
cNoseSurface = gmsh.model.geo.addPlaneSurface([cNoseLoop])
cTopSurface = gmsh.model.geo.addPlaneSurface([cTopLoop])
cBottomSurface = gmsh.model.geo.addPlaneSurface([cBottomLoop])
gmsh.model.geo.mesh.setTransfiniteSurface(cNoseSurface, "Left", [cBottomLeftPoint, teLowerPoints[0], lePoints[0], cTopLeftPoint])
gmsh.model.geo.mesh.setTransfiniteSurface(cTopSurface, "Left", [cTopRightPoint, cTopLeftPoint, lePoints[0], cRightPoint])
gmsh.model.geo.mesh.setRecombine(2, cNoseSurface)
gmsh.model.geo.mesh.setRecombine(2, cTopSurface)
gmsh.model.geo.mesh.setRecombine(2, cBottomSurface)

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