import gmsh
import sys
import numpy as np
import math

################################################################################
#   Input Parameters
################################################################################

aerofoilName = "NACA0009"
Re = 80000
alpha = 15

ffBoundaryScaling = 40
blNumPoints = 75
cNoseNearWallHeight = 0.001
outerLayerHeight = 1
outerLayerLength = 11
outerLayerCellSizeScaling = 3

blThicknessProgression = 1.2
blTEProgression = 1.04
blLEBump = 10
cProgression = 1.2
wakeProgression = 1.2
blPointsSplit = 0.45
cMeshHeight = 0.2
cWakeLength = 1

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
deltaLam = 1.227/(Re**0.75)
firstCellHeight = deltaLam/0.047
print("First cell height laminar flow: "+str(firstCellHeight))

cNoseNearWallHeight = firstCellHeight*0.2

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
#   C Mesh Geometry
################################################################################

#Create C mesh points
cTopLeftPoint = gmsh.model.geo.addPoint(thickestXCoord, cMeshHeight, 0)
cLeftPoint = gmsh.model.geo.addPoint(-cNoseDistance, 0, 0)
cBottomLeftPoint = gmsh.model.geo.addPoint(thickestXCoord, -cMeshHeight, 0)
cTopRightPoint = gmsh.model.geo.addPoint(wakeProgressionLength+1, cMeshHeight, 0)
cRightPoint = gmsh.model.geo.addPoint(wakeProgressionLength+1, 0, 0)
cBottomRightPoint = gmsh.model.geo.addPoint(wakeProgressionLength+1, -cMeshHeight, 0)
cFarRightTopPoint  = gmsh.model.geo.addPoint(cWakeLength+1, cMeshHeight, 0)
cFarRightBottomPoint  = gmsh.model.geo.addPoint(cWakeLength+1, -cMeshHeight, 0)

#Create C mesh lines
cTopArc = gmsh.model.geo.addEllipseArc(cTopLeftPoint, thickestX, cLeftPoint, cLeftPoint)
cBottomArc = gmsh.model.geo.addEllipseArc(cLeftPoint, thickestX, cLeftPoint, cBottomLeftPoint)
cTopLine = gmsh.model.geo.addLine(cTopRightPoint, cTopLeftPoint)
cBottomLine = gmsh.model.geo.addLine(cBottomLeftPoint, cBottomRightPoint)
cTopRightLine = gmsh.model.geo.addLine(cRightPoint, cTopRightPoint)
cBottomRightLine = gmsh.model.geo.addLine(cBottomRightPoint, cRightPoint)
cMidLine = gmsh.model.geo.addLine(trailingEdge, cRightPoint)
cFarRightTopLine = gmsh.model.geo.addLine(cFarRightTopPoint, cTopRightPoint)
cFarRightBottomLine = gmsh.model.geo.addLine(cBottomRightPoint, cFarRightBottomPoint)
cFarRightLine = gmsh.model.geo.addLine(cFarRightBottomPoint, cFarRightTopPoint)

#Set transfinite C mesh lines
gmsh.model.geo.mesh.setTransfiniteCurve(cTopArc, blLENumPoints+1)
gmsh.model.geo.mesh.setTransfiniteCurve(cBottomArc, blLENumPoints+1)
gmsh.model.geo.mesh.setTransfiniteCurve(cMidLine, wakeNumPoints, "Progression", wakeProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(cTopLine, wakeNumPoints+blTENumPoints-1)
gmsh.model.geo.mesh.setTransfiniteCurve(cBottomLine, wakeNumPoints+blTENumPoints-1)
gmsh.model.geo.mesh.setTransfiniteCurve(cTopRightLine, cNumPoints, "Progression", cProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(cBottomRightLine, cNumPoints, "Progression", -cProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(cFarRightLine, cNumPoints*2-1)
gmsh.model.geo.mesh.setTransfiniteCurve(cFarRightTopLine, round((cWakeLength-wakeProgressionLength)/cFarCellSize))
gmsh.model.geo.mesh.setTransfiniteCurve(cFarRightBottomLine, round((cWakeLength-wakeProgressionLength)/cFarCellSize))

#Create C mesh curve loop
cLoop = gmsh.model.geo.addCurveLoop([cTopLine, cTopArc, cBottomArc, cBottomLine, cBottomRightLine, -cMidLine, -teLowerSpline, -leSpline, -teUpperSpline, cMidLine, cTopRightLine])
cFarLoop = gmsh.model.geo.addCurveLoop([cFarRightTopLine, -cTopRightLine, -cBottomRightLine, cFarRightBottomLine, cFarRightLine])

#Create C mesh surface
cSurface = gmsh.model.geo.addPlaneSurface([cLoop])
cFarSurface = gmsh.model.geo.addPlaneSurface([cFarLoop])
gmsh.model.geo.mesh.setTransfiniteSurface(cSurface, "Left", [cRightPoint, cTopRightPoint, cBottomRightPoint, cRightPoint])
gmsh.model.geo.mesh.setTransfiniteSurface(cFarSurface, "Left", [cFarRightTopPoint, cTopRightPoint, cBottomRightPoint, cFarRightBottomPoint])
gmsh.model.geo.mesh.setRecombine(2, cSurface)
gmsh.model.geo.mesh.setRecombine(2, cFarSurface)

#Rotate C mesh for desired alpha
gmsh.model.geo.rotate([(2, cSurface),(2, cFarSurface)], 1, 0, 0, 0, 0, 1, math.radians(-alpha))

################################################################################
#   Outer Layer Geometry
################################################################################

#Create outer layer mesh points
olTopLeftPoint = gmsh.model.geo.addPoint(0, outerLayerHeight, 0, cFarCellSize*outerLayerCellSizeScaling)
olLeftPoint = gmsh.model.geo.addPoint(-outerLayerHeight, 0, 0, cFarCellSize*outerLayerCellSizeScaling)
olBottomLeftPoint = gmsh.model.geo.addPoint(0, -outerLayerHeight, 0, cFarCellSize*outerLayerCellSizeScaling)
olTopRightPoint = gmsh.model.geo.addPoint(outerLayerLength, outerLayerHeight, 0, cFarCellSize*outerLayerCellSizeScaling)
olBottomRightPoint = gmsh.model.geo.addPoint(outerLayerLength, -outerLayerHeight, 0, cFarCellSize*outerLayerCellSizeScaling)

#Create outer layer mesh lines
olTopArc = gmsh.model.geo.addCircleArc(olTopLeftPoint, origin, olLeftPoint)
olBottomArc = gmsh.model.geo.addCircleArc(olLeftPoint, origin, olBottomLeftPoint)
olTopLine = gmsh.model.geo.addLine(olTopRightPoint, olTopLeftPoint)
olBottomLine = gmsh.model.geo.addLine(olBottomLeftPoint, olBottomRightPoint)
olRightLine = gmsh.model.geo.addLine(olBottomRightPoint, olTopRightPoint)

#Create outer layer mesh curve loop
olLoop = gmsh.model.geo.addCurveLoop([olTopLine, olTopArc, olBottomArc, olBottomLine, olRightLine])

#Create outer layer mesh surface
olSurface = gmsh.model.geo.addPlaneSurface([olLoop, cLoop, cFarLoop])
#gmsh.model.geo.mesh.setRecombine(2, olSurface)

################################################################################
#   Far Field Geometry
################################################################################

#Create points for the far field
ffTopPoint = gmsh.model.geo.addPoint(1, ffBoundaryScaling, 0)
ffLeftPoint = gmsh.model.geo.addPoint(1-ffBoundaryScaling, 0, 0)
ffBottomPoint = gmsh.model.geo.addPoint(1, -ffBoundaryScaling, 0)
ffRightPoint = gmsh.model.geo.addPoint(1+ffBoundaryScaling, 0, 0)

#Create lines for the far field
ffTopLeftArc = gmsh.model.geo.addCircleArc(ffTopPoint, trailingEdge, ffLeftPoint)
ffBottomLeftArc = gmsh.model.geo.addCircleArc(ffLeftPoint, trailingEdge, ffBottomPoint)
ffTopRightArc = gmsh.model.geo.addCircleArc(ffBottomPoint, trailingEdge, ffRightPoint)
ffBottomRightArc = gmsh.model.geo.addCircleArc(ffRightPoint, trailingEdge, ffTopPoint)

#Create the curve loop for the far field
farFieldLoop = gmsh.model.geo.addCurveLoop([ffTopLeftArc, ffBottomLeftArc, ffBottomRightArc, ffTopRightArc])

#Create the surface
farFieldSurface = gmsh.model.geo.addPlaneSurface([farFieldLoop, olSurface])

################################################################################
#   Synchronize
################################################################################

#Synchronize CAD entities with the Gmsh model
gmsh.model.geo.synchronize()

################################################################################
#   Physical Groups
################################################################################

#Create physical groups
aerofoilGroup = gmsh.model.addPhysicalGroup(1, [teUpperSpline, leSpline, teLowerSpline])
gmsh.model.setPhysicalName(1, aerofoilGroup, "aerofoilwall")
farFieldGroup = gmsh.model.addPhysicalGroup(1, [ffTopLeftArc, ffBottomLeftArc, ffBottomRightArc, ffTopRightArc])
gmsh.model.setPhysicalName(1, farFieldGroup, "farfieldboundary")
fluidGroup = gmsh.model.addPhysicalGroup(2, [cSurface, cFarSurface, olSurface, farFieldSurface])
gmsh.model.setPhysicalName(2, fluidGroup, "fluid")

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