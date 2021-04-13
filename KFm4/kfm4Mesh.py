import gmsh
import sys
import numpy as np
import math

################################################################################
#   Input Parameters
################################################################################

aerofoilName = "NACA0009"
Re = 80000
alpha = 20

ffBoundaryScaling = 80
leNumPoints = 25
teNumPoints = 25
upperNumPoints = 12
cNoseNearWallHeight = 0.001
outerLayerHeight = 1
outerLayerLength = 11
outerLayerCellSizeScaling = 3

blThicknessProgression = 1.2
blLEBump = 10
cProgression = 1.2
wakeProgression = 1.2
blPointsSplit = 0.45
cMeshHeight = 0.2
cWakeLength = 1

ffDensity = 10

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

#Calculate the the first cell height required to resolve laminar flow 
deltaLam = 1.3016/(Re**0.75)
firstCellHeight = deltaLam/0.047
print("First cell height laminar flow: "+str(firstCellHeight))

cNoseNearWallHeight = firstCellHeight*0.25

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
teCellWidth = firstCellHeight
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
for i in range(0, len(coords)-1):
    aerofoilPoints.append(gmsh.model.geo.addPoint(coords[i][0], coords[i][1], 0))
aerofoilPoints.append(aerofoilPoints[0])
upperStepPoint = gmsh.model.geo.addPoint(0.5, coords[24][1], 0)
lowerStepPoint = gmsh.model.geo.addPoint(0.5, -coords[24][1], 0)

#Divide up the aerofoil point into leading and upper and lower trailing edges
teUpperPoints = []
upperPoints = []
lePoints = []
lowerPoints = []
teLowerPoints = []
for i in range(0, 25):
    teUpperPoints.append(aerofoilPoints[i])
for i in range(50, thickestIndex+1):
    upperPoints.append(aerofoilPoints[i])
for i in range(thickestIndex, 201-thickestIndex):
    lePoints.append(aerofoilPoints[i])
for i in range(176, 201):
    teLowerPoints.append(aerofoilPoints[i])
for i in range(200-thickestIndex, 151):
    lowerPoints.append(aerofoilPoints[i])

#Create the aerofoil lines
teUpperSpline = gmsh.model.geo.addSpline(teUpperPoints)
leSpline = gmsh.model.geo.addSpline(lePoints)
upperSpline = gmsh.model.geo.addSpline(upperPoints)
lowerSpline = gmsh.model.geo.addSpline(lowerPoints)
teLowerSpline = gmsh.model.geo.addSpline(teLowerPoints)
stepUpperHorzLine = gmsh.model.geo.addLine(aerofoilPoints[24], upperStepPoint)
stepLowerHorzLine = gmsh.model.geo.addLine(lowerStepPoint, aerofoilPoints[176])
stepUpperVertLine = gmsh.model.geo.addLine(upperStepPoint, aerofoilPoints[50])
stepLowerVertLine = gmsh.model.geo.addLine(aerofoilPoints[150], lowerStepPoint)

#Set aerofoil transfinite curves
gmsh.model.geo.mesh.setTransfiniteCurve(leSpline, 2*leNumPoints-1, "Bump", blLEBump)
gmsh.model.geo.mesh.setTransfiniteCurve(stepUpperHorzLine, int(round(0.4*(coords[24][0]-0.5)/firstCellHeight)), "Bump", 0.4)
gmsh.model.geo.mesh.setTransfiniteCurve(stepLowerHorzLine, int(round(0.4*(coords[24][0]-0.5)/firstCellHeight)), "Bump", 0.4)
gmsh.model.geo.mesh.setTransfiniteCurve(stepUpperVertLine, int(round((coords[50][1])/firstCellHeight)))
gmsh.model.geo.mesh.setTransfiniteCurve(stepLowerVertLine, int(round((coords[50][1])/firstCellHeight)))
gmsh.model.geo.mesh.setTransfiniteCurve(teUpperSpline, int(round(0.5*(1-coords[24][0])/firstCellHeight)), "Bump", 0.4)
gmsh.model.geo.mesh.setTransfiniteCurve(teLowerSpline, int(round(0.5*(1-coords[24][0])/firstCellHeight)), "Bump", 0.4)
gmsh.model.geo.mesh.setTransfiniteCurve(upperSpline, upperNumPoints, "Progression", 1.15)
gmsh.model.geo.mesh.setTransfiniteCurve(lowerSpline, upperNumPoints, "Progression", -1.15)

################################################################################
#   Step Mesh Geometry
################################################################################

#Create step mesh points
stepTopPoint = gmsh.model.geo.addPoint(coords[24][0], coords[50][1], 0)
stepTopRightPoint = gmsh.model.geo.addPoint(1, coords[50][1], 0)
stepTopFarRightPoint = gmsh.model.geo.addPoint(wakeProgressionLength+1, coords[50][1], 0)
stepBottomPoint = gmsh.model.geo.addPoint(coords[24][0], -coords[50][1], 0)
stepBottomRightPoint = gmsh.model.geo.addPoint(1, -coords[50][1], 0)
stepBottomFarRightPoint = gmsh.model.geo.addPoint(wakeProgressionLength+1, -coords[50][1], 0)

#Create step mesh lines
stepTopLine = gmsh.model.geo.addLine(stepTopPoint, aerofoilPoints[50])
stepUpperMidVertLine = gmsh.model.geo.addLine(aerofoilPoints[24], stepTopPoint)
stepBottomLine = gmsh.model.geo.addLine(aerofoilPoints[150], stepBottomPoint)
stepLowerMidVertLine = gmsh.model.geo.addLine(stepBottomPoint, aerofoilPoints[176])
stepTopRightLine = gmsh.model.geo.addLine(stepTopRightPoint, stepTopPoint)
stepBottomRightLine = gmsh.model.geo.addLine(stepBottomPoint, stepBottomRightPoint)
stepUpperMidRightVertLine = gmsh.model.geo.addLine(aerofoilPoints[0], stepTopRightPoint)
stepLowerMidRightVertLine = gmsh.model.geo.addLine(stepBottomRightPoint, aerofoilPoints[-1])
stepFarRightLine = gmsh.model.geo.addLine(stepBottomFarRightPoint, stepTopFarRightPoint)
stepUpperFarRightLine = gmsh.model.geo.addLine(stepTopFarRightPoint, stepTopRightPoint)
stepLowerFarRightLine = gmsh.model.geo.addLine(stepBottomRightPoint, stepBottomFarRightPoint)

#Set transfinite step mesh lines
gmsh.model.geo.mesh.setTransfiniteCurve(stepTopLine, int(round(0.4*(coords[24][0]-0.5)/firstCellHeight)), "Bump", 0.4)
gmsh.model.geo.mesh.setTransfiniteCurve(stepBottomLine, int(round(0.4*(coords[24][0]-0.5)/firstCellHeight)), "Bump", 0.4)
gmsh.model.geo.mesh.setTransfiniteCurve(stepUpperMidVertLine, int(round((coords[50][1])/firstCellHeight)))
gmsh.model.geo.mesh.setTransfiniteCurve(stepLowerMidVertLine, int(round((coords[50][1])/firstCellHeight)))
gmsh.model.geo.mesh.setTransfiniteCurve(stepUpperMidRightVertLine, int(round((coords[50][1])/firstCellHeight)))
gmsh.model.geo.mesh.setTransfiniteCurve(stepLowerMidRightVertLine, int(round((coords[50][1])/firstCellHeight)))
gmsh.model.geo.mesh.setTransfiniteCurve(stepTopRightLine, int(round(0.5*(1-coords[24][0])/firstCellHeight)), "Bump", 0.4)
gmsh.model.geo.mesh.setTransfiniteCurve(stepBottomRightLine, int(round(0.5*(1-coords[24][0])/firstCellHeight)), "Bump", 0.4)
gmsh.model.geo.mesh.setTransfiniteCurve(stepFarRightLine, 2*int(round((coords[50][1])/firstCellHeight))-1)
gmsh.model.geo.mesh.setTransfiniteCurve(stepUpperFarRightLine, wakeNumPoints, "Progression", -wakeProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(stepLowerFarRightLine, wakeNumPoints, "Progression", wakeProgression)

#Create step mesh curve loops
stepTopLoop = gmsh.model.geo.addCurveLoop([stepTopLine, -stepUpperVertLine, -stepUpperHorzLine, stepUpperMidVertLine])
stepTopRightLoop = gmsh.model.geo.addCurveLoop([stepTopRightLine, -stepUpperMidVertLine, -teUpperSpline, stepUpperMidRightVertLine])
stepBottomLoop = gmsh.model.geo.addCurveLoop([-stepLowerHorzLine, -stepLowerVertLine, stepBottomLine, stepLowerMidVertLine])
stepBottomRightLoop = gmsh.model.geo.addCurveLoop([-teLowerSpline, -stepLowerMidVertLine, stepBottomRightLine, stepLowerMidRightVertLine])
stepFarRightLoop = gmsh.model.geo.addCurveLoop([stepUpperFarRightLine, -stepUpperMidRightVertLine, -stepLowerMidRightVertLine, stepLowerFarRightLine, stepFarRightLine])

#Create step mesh surface
stepTopSurface = gmsh.model.geo.addPlaneSurface([stepTopLoop])
stepTopRightSurface = gmsh.model.geo.addPlaneSurface([stepTopRightLoop])
stepBottomSurface = gmsh.model.geo.addPlaneSurface([stepBottomLoop])
stepBottomRightSurface = gmsh.model.geo.addPlaneSurface([stepBottomRightLoop])
stepFarRightSurface = gmsh.model.geo.addPlaneSurface([stepFarRightLoop])
gmsh.model.geo.mesh.setTransfiniteSurface(stepTopSurface)
gmsh.model.geo.mesh.setTransfiniteSurface(stepBottomSurface)
gmsh.model.geo.mesh.setTransfiniteSurface(stepTopRightSurface)
gmsh.model.geo.mesh.setTransfiniteSurface(stepBottomRightSurface)
gmsh.model.geo.mesh.setTransfiniteSurface(stepFarRightSurface, "Left", [stepTopFarRightPoint, stepTopRightPoint, stepBottomRightPoint, stepBottomFarRightPoint])
gmsh.model.geo.mesh.setRecombine(2, stepTopSurface)
gmsh.model.geo.mesh.setRecombine(2, stepBottomSurface)
gmsh.model.geo.mesh.setRecombine(2, stepTopRightSurface)
gmsh.model.geo.mesh.setRecombine(2, stepBottomRightSurface)
gmsh.model.geo.mesh.setRecombine(2, stepFarRightSurface)

################################################################################
#   C Mesh Geometry
################################################################################

#Create C mesh points
cTopLeftPoint = gmsh.model.geo.addPoint(thickestXCoord, cMeshHeight, 0)
cLeftPoint = gmsh.model.geo.addPoint(-cNoseDistance, 0, 0)
cBottomLeftPoint = gmsh.model.geo.addPoint(thickestXCoord, -cMeshHeight, 0)
cTopRightPoint = gmsh.model.geo.addPoint(wakeProgressionLength+1, cMeshHeight, 0)
cBottomRightPoint = gmsh.model.geo.addPoint(wakeProgressionLength+1, -cMeshHeight, 0)
cFarRightTopPoint  = gmsh.model.geo.addPoint(cWakeLength+1, cMeshHeight, 0)
cFarRightBottomPoint  = gmsh.model.geo.addPoint(cWakeLength+1, -cMeshHeight, 0)

#Create C mesh lines
cTopArc = gmsh.model.geo.addEllipseArc(cTopLeftPoint, thickestX, cLeftPoint, cLeftPoint)
cBottomArc = gmsh.model.geo.addEllipseArc(cLeftPoint, thickestX, cLeftPoint, cBottomLeftPoint)
cTopLine = gmsh.model.geo.addLine(cTopRightPoint, cTopLeftPoint)
cBottomLine = gmsh.model.geo.addLine(cBottomLeftPoint, cBottomRightPoint)
cTopRightLine = gmsh.model.geo.addLine(stepTopFarRightPoint, cTopRightPoint)
cBottomRightLine = gmsh.model.geo.addLine(cBottomRightPoint, stepBottomFarRightPoint)
cFarRightTopLine = gmsh.model.geo.addLine(cFarRightTopPoint, cTopRightPoint)
cFarRightBottomLine = gmsh.model.geo.addLine(cBottomRightPoint, cFarRightBottomPoint)
cFarRightLine = gmsh.model.geo.addLine(cFarRightBottomPoint, cFarRightTopPoint)

#Set transfinite C mesh lines
gmsh.model.geo.mesh.setTransfiniteCurve(cTopArc, leNumPoints)
gmsh.model.geo.mesh.setTransfiniteCurve(cBottomArc, leNumPoints)
gmsh.model.geo.mesh.setTransfiniteCurve(cTopLine, upperNumPoints+wakeNumPoints+int(round(0.4*(coords[24][0]-0.5)/firstCellHeight))+int(round(0.5*(1-coords[24][0])/firstCellHeight))-3)
gmsh.model.geo.mesh.setTransfiniteCurve(cBottomLine, upperNumPoints+wakeNumPoints+int(round(0.4*(coords[24][0]-0.5)/firstCellHeight))+int(round(0.5*(1-coords[24][0])/firstCellHeight))-3)
gmsh.model.geo.mesh.setTransfiniteCurve(cTopRightLine, cNumPoints, "Progression", cProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(cBottomRightLine, cNumPoints, "Progression", -cProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(cFarRightLine, 2*int(round((coords[50][1])/firstCellHeight))+cNumPoints*2-3)
gmsh.model.geo.mesh.setTransfiniteCurve(cFarRightTopLine, int(round((cWakeLength-wakeProgressionLength)/cFarCellSize)))
gmsh.model.geo.mesh.setTransfiniteCurve(cFarRightBottomLine, int(round((cWakeLength-wakeProgressionLength)/cFarCellSize)))

#Create C mesh curve loop
cLoop = gmsh.model.geo.addCurveLoop([cTopLine, cTopArc, cBottomArc, cBottomLine, cBottomRightLine, -stepLowerFarRightLine, -stepBottomRightLine, -stepBottomLine, -lowerSpline, -leSpline, -upperSpline, -stepTopLine, -stepTopRightLine, -stepUpperFarRightLine, cTopRightLine])
cFarLoop = gmsh.model.geo.addCurveLoop([cFarRightTopLine, -cTopRightLine, -stepFarRightLine, -cBottomRightLine, cFarRightBottomLine, cFarRightLine])

#Create C mesh surface
cSurface = gmsh.model.geo.addPlaneSurface([cLoop])
cFarSurface = gmsh.model.geo.addPlaneSurface([cFarLoop])
gmsh.model.geo.mesh.setTransfiniteSurface(cSurface, "Left", [stepTopFarRightPoint, cTopRightPoint, cBottomRightPoint, stepBottomFarRightPoint])
gmsh.model.geo.mesh.setTransfiniteSurface(cFarSurface, "Left", [cFarRightTopPoint, cTopRightPoint, cBottomRightPoint, cFarRightBottomPoint])
gmsh.model.geo.mesh.setRecombine(2, cSurface)
gmsh.model.geo.mesh.setRecombine(2, cFarSurface)

#Rotate C mesh for desired alpha
gmsh.model.geo.rotate([(2, cSurface),(2, cFarSurface), (2, stepTopSurface), (2, stepBottomSurface), (2, stepTopRightSurface), (2, stepBottomRightSurface)], 1, 0, 0, 0, 0, 1, math.radians(-alpha))

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
ffTopPoint = gmsh.model.geo.addPoint(1, ffBoundaryScaling, 0, ffDensity)
ffLeftPoint = gmsh.model.geo.addPoint(1-ffBoundaryScaling, 0, 0, ffDensity)
ffBottomPoint = gmsh.model.geo.addPoint(1, -ffBoundaryScaling, 0, ffDensity)
ffRightPoint = gmsh.model.geo.addPoint(1+ffBoundaryScaling, 0, 0, ffDensity)

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
fluidGroup = gmsh.model.addPhysicalGroup(2, [cSurface, cFarSurface, olSurface, farFieldSurface, stepTopSurface,stepBottomSurface, stepTopRightSurface, stepBottomRightSurface, stepFarRightSurface])
gmsh.model.setPhysicalName(2, fluidGroup, "fluid")

################################################################################
#   Output
################################################################################

#Set colours
# for i in range(0, len(boundaryLayerQuadSurfaces)):
#     gmsh.model.setColor([(2, boundaryLayerQuadSurfaces[i])], 255, 0, 0)
# gmsh.model.setColor([(2,volumeSurface)], 0, 255, 0)
gmsh.option.setNumber("General.BackgroundGradient", 0)

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