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
blThicknessNumPoints = 3
aerofoilThickness = 0.09
Re = 100000
cMeshHeight = 1

srAngle = 30
srThickness = 0.4
blThicknessProgression = 1.2
blTEProgression = 1.07
blLEBump = 6
ffProgression = 1.2
wakeAngle = 2
wakeProgression = 1.08

################################################################################
#   Calculate Useful Parameters
################################################################################

#Boundary layer parameters
blThickness = 5/math.sqrt(Re)
print("Theoretical boundary layer thickness: "+str(blThickness))

deltaLam = 1.227/(Re**0.75)
firstCellHeightLam = deltaLam/0.0695
print("First cell heigh laminar flow: "+str(firstCellHeightLam))

deltaTurb = (13.1463**0.875)/(Re**0.9)
firstCellHeightTurb = deltaTurb/0.0695
print("First cell heigh turbulent flow: "+str(firstCellHeightTurb))

nearWallCellThickness = firstCellHeightLam
print("Cell size near wall: "+str(nearWallCellThickness))

blThickness = nearWallCellThickness*(1-blThicknessProgression**blThicknessNumPoints)/(1-blThicknessProgression)
print("Mesh boundary layer thickness: "+str(blThickness))

#Wake parameters
wakeLength = 2*ffBoundaryScaling
wakeEndThickness = (ffBoundaryScaling*2-1)*math.tan(math.radians(wakeAngle))

teCellWidth = 0.7*(1-blTEProgression)/(1-blTEProgression**(blNumPoints*0.7-1))
NWake = round(math.log(1-(wakeLength*(1-wakeProgression)/teCellWidth))/math.log(wakeProgression))

#Stall region number of points
srFirstThickness = nearWallCellThickness*blThicknessProgression**(blThicknessNumPoints-1)
srProgression = blThicknessProgression
srNumPoints = round(math.log(1-((srThickness-blThickness-aerofoilThickness/2)*(1-srProgression)/srFirstThickness))/math.log(srProgression))

#Density of far field points
ffPointsPerLength = 1/(nearWallCellThickness*ffBoundaryScaling*40)
ffPointsPerLength = 5/ffBoundaryScaling

################################################################################
#   Initialize Gmsh
################################################################################

#Initialize Gmsh and name the model
gmsh.initialize()
gmsh.model.add(aerofoilName)

################################################################################
#   Aerofoil Geometry
################################################################################

#Get the coordinates of the aerofoil
aerofoilData = np.loadtxt(aerofoilName+".dat")
coords = aerofoilData

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
blTENumPoints = round((1-0.35)*blNumPoints)
blLENumPoints = round(0.35*blNumPoints)*2+1
teUpperSpline = gmsh.model.geo.addSpline(teUpperPoints)
leSpline = gmsh.model.geo.addSpline(lePoints)
teLowerSpline = gmsh.model.geo.addSpline(teLowerPoints)
gmsh.model.geo.mesh.setTransfiniteCurve(teUpperSpline, blTENumPoints, "Progression", blTEProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(leSpline, blLENumPoints, "Bump", blLEBump)
gmsh.model.geo.mesh.setTransfiniteCurve(teLowerSpline, blTENumPoints, "Progression", -blTEProgression)

################################################################################
#   C Mesh Geometry
################################################################################

#Create C mesh points
cNosePoints = []
for i in range(thickestIndex, len(coords)-thickestIndex):
    cNosePoints.append(gmsh.model.geo.addPoint(coords[i][0], coords[i][1], 0))

#Create C mesh lines
cSpline = gmsh.model.geo.addSpline(cNosePoints)
gmsh.model.geo.dilate([(1, cSpline)], coords[thickestIndex][0], 0, 0, 1+(cMeshHeight)/coords[thickestIndex][0], (coords[thickestIndex][1]+cMeshHeight)/coords[thickestIndex][1], 1)


################################################################################
#   Boundary Layer Geometry
################################################################################

# #Create boundary layer points
# blTopRightPoint = gmsh.model.geo.addPoint(1, blThickness, 0)
# blBottomRightPoint = gmsh.model.geo.addPoint(1, -blThickness, 0)
# blTopMidPoint = gmsh.model.geo.addPoint(coords[thickestIndex][0], coords[thickestIndex][1]+blThickness, 0)
# blBottomMidPoint = gmsh.model.geo.addPoint(coords[thickestIndex][0], -coords[thickestIndex][1]-blThickness, 0)
# blNosePoints = []
# blTeUpperPoints = []
# blTeLowerPoints = []
# for i in range(thickestIndex, len(coords)-thickestIndex):
#     blNosePoints.append(gmsh.model.geo.addPoint(coords[i][0], coords[i][1], 0))
# for i in range(0, thickestIndex+1):
#     blTeUpperPoints.append(gmsh.model.geo.addPoint(coords[i][0], coords[i][1]+blThickness, 0))
# for i in range(len(coords)-thickestIndex-1, len(coords)):
#     blTeLowerPoints.append(gmsh.model.geo.addPoint(coords[i][0], coords[i][1]-blThickness, 0))

# #Create the boundary layer lines
# #blTopTELine = gmsh.model.geo.addLine(blTopRightPoint, blTopMidPoint)
# #blBottomTELine = gmsh.model.geo.addLine(blBottomMidPoint, blBottomRightPoint)
# blTopRightLine = gmsh.model.geo.addLine(aerofoilPoints[0], blTopRightPoint)
# blBottomRightLine = gmsh.model.geo.addLine(blBottomRightPoint, aerofoilPoints[-1])
# blTopMidLine = gmsh.model.geo.addLine(blTopMidPoint, aerofoilPoints[thickestIndex])
# blBottomMidLine = gmsh.model.geo.addLine(aerofoilPoints[len(coords)-thickestIndex-1], blBottomMidPoint)
# blLESpline = gmsh.model.geo.addSpline(blNosePoints)
# blTEUpperSpline = gmsh.model.geo.addSpline(blTeUpperPoints)
# blTELowerSpline = gmsh.model.geo.addSpline(blTeLowerPoints)
# gmsh.model.geo.dilate([(1, blLESpline)], coords[thickestIndex][0], 0, 0, 1+(2*blThickness)/coords[thickestIndex][0], (coords[thickestIndex][1]+blThickness)/coords[thickestIndex][1], 1)
# gmsh.model.geo.mesh.setTransfiniteCurve(blTEUpperSpline, blTENumPoints, "Progression", blTEProgression)
# gmsh.model.geo.mesh.setTransfiniteCurve(blTELowerSpline, blTENumPoints, "Progression", -blTEProgression)
# gmsh.model.geo.mesh.setTransfiniteCurve(blTopRightLine, blThicknessNumPoints, "Progression", blThicknessProgression)
# gmsh.model.geo.mesh.setTransfiniteCurve(blBottomRightLine, blThicknessNumPoints, "Progression", -blThicknessProgression)
# gmsh.model.geo.mesh.setTransfiniteCurve(blTopMidLine, blThicknessNumPoints, "Progression", -blThicknessProgression)
# gmsh.model.geo.mesh.setTransfiniteCurve(blBottomMidLine, blThicknessNumPoints, "Progression", blThicknessProgression)
# gmsh.model.geo.mesh.setTransfiniteCurve(blLESpline, blLENumPoints, "Bump", blLEBump)

# #Create boundary layer curve loops
# blTopTELoop = gmsh.model.geo.addCurveLoop([blTEUpperSpline, blTopMidLine, -teUpperSpline, blTopRightLine])
# blLELoop = gmsh.model.geo.addCurveLoop([blLESpline, -blBottomMidLine, -leSpline, -blTopMidLine])
# blBottomTELoop = gmsh.model.geo.addCurveLoop([-teLowerSpline, blBottomMidLine, blTELowerSpline, blBottomRightLine])

# #Create boundary layer surfaces
# blTopTESurface = gmsh.model.geo.addPlaneSurface([blTopTELoop])
# blLESurface = gmsh.model.geo.addPlaneSurface([blLELoop])
# blBottomTESurface = gmsh.model.geo.addPlaneSurface([blBottomTELoop])
# gmsh.model.geo.mesh.setTransfiniteSurface(blTopTESurface)
# gmsh.model.geo.mesh.setTransfiniteSurface(blLESurface)
# gmsh.model.geo.mesh.setTransfiniteSurface(blBottomTESurface)
# gmsh.model.geo.mesh.setRecombine(2, blTopTESurface)
# gmsh.model.geo.mesh.setRecombine(2, blLESurface)
# gmsh.model.geo.mesh.setRecombine(2, blBottomTESurface)

################################################################################
#   Wake Geometry
################################################################################

# #Create points for wake
# wakeTopLeftPoint = blTopRightPoint
# wakeBottomLeftPoint = blBottomRightPoint
# wakeTopRightPoint = gmsh.model.geo.addPoint(wakeLength, wakeEndThickness, 0)
# wakeBottomRightPoint = gmsh.model.geo.addPoint(wakeLength, -wakeEndThickness, 0)

# #Create lines for wake
# wakeTopLine = gmsh.model.geo.addLine(wakeTopRightPoint, wakeTopLeftPoint)
# wakeBottomLine = gmsh.model.geo.addLine(wakeBottomLeftPoint, wakeBottomRightPoint)
# wakeRightLine = gmsh.model.geo.addLine(wakeBottomRightPoint, wakeTopRightPoint)
# gmsh.model.geo.mesh.setTransfiniteCurve(wakeTopLine, NWake, "Progression", -wakeProgression)
# gmsh.model.geo.mesh.setTransfiniteCurve(wakeBottomLine, NWake, "Progression", wakeProgression)
# gmsh.model.geo.mesh.setTransfiniteCurve(wakeRightLine, blThicknessNumPoints*2-1)

# #Create curveloop for wake
# wakeLoop = gmsh.model.geo.addCurveLoop([wakeTopLine, -blTopRightLine, -blBottomRightLine, wakeBottomLine, wakeRightLine])

# #Create surface for wake
# wakeSurface = gmsh.model.geo.addPlaneSurface([wakeLoop])
# gmsh.model.geo.mesh.setTransfiniteSurface(wakeSurface, "Left", [wakeTopRightPoint, wakeTopLeftPoint, wakeBottomLeftPoint, wakeBottomRightPoint])
# gmsh.model.geo.mesh.setRecombine(2, wakeSurface)


################################################################################
#   Far Field Geometry
################################################################################

# #Create points for the far field
# ffTopRightPoint = gmsh.model.geo.addPoint(wakeLength, ffBoundaryScaling, 0)
# ffTopLeftPoint = gmsh.model.geo.addPoint(-ffBoundaryScaling, ffBoundaryScaling, 0)
# ffBottomLeftPoint = gmsh.model.geo.addPoint(-ffBoundaryScaling, -ffBoundaryScaling, 0)
# ffBottomRightPoint = gmsh.model.geo.addPoint(wakeLength, -ffBoundaryScaling, 0)

# #Create lines for the far field
# ffLeftLine = gmsh.model.geo.addLine(ffTopLeftPoint, ffBottomLeftPoint)
# ffBottomLine = gmsh.model.geo.addLine(ffBottomLeftPoint, ffBottomRightPoint)
# ffBottomRightLine = gmsh.model.geo.addLine(ffBottomRightPoint, wakeBottomRightPoint)
# ffRightLine2 = gmsh.model.geo.addLine(wakeTopRightPoint, ffTopRightPoint)
# gmsh.model.geo.mesh.setTransfiniteCurve(ffLeftLine, round(ffPointsPerLength*ffBoundaryScaling*2))
# gmsh.model.geo.mesh.setTransfiniteCurve(ffBottomLine, round(ffPointsPerLength*(ffBoundaryScaling+wakeLength)))
# gmsh.model.geo.mesh.setTransfiniteCurve(ffBottomRightLine, round(ffPointsPerLength*(ffBoundaryScaling-wakeEndThickness/2)))
# gmsh.model.geo.mesh.setTransfiniteCurve(ffRightLine2, round(ffPointsPerLength*(ffBoundaryScaling-wakeEndThickness/2)))

# #Create the curve loop for the far field
# farFieldLoop = gmsh.model.geo.addCurveLoop([ffTopLine, ffLeftLine, ffBottomLine, ffBottomRightLine, -wakeBottomLine, -blTELowerSpline, -srBottomLine, -srSpline, -srMidLine, -srTopLine])

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