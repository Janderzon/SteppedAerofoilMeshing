import gmsh
import sys
import numpy as np
import math

################################################################################
#   Input Parameters
################################################################################

blTEProgression = 1.1
blTENumPoints = 10
blLEBump = 5
blLENumPoints = 21
blThickness = 0.2
blThicknessProgression = 1.2
wakeLength = 20
wakeEndThickness = 3
wakeNumPoints = 40
wakeProgression = 1.1
ffHeight = 20
ffHorzOffset = 10
ffPointsPerLength = 0.5
aerofoilName = "NACA63-009"
nearWallCellThickness = 0.02

################################################################################
#   Calculate Useful Parameters
################################################################################

N = math.log(1-blThickness*(1-blThicknessProgression)/nearWallCellThickness)/math.log(blThicknessProgression)
blThicknessNumPoints = round(N)+1
nearWallCellThicknessRounded = blThickness*(1-blThicknessProgression)/(1-blThicknessProgression**blThicknessNumPoints)
print("Cell size near wall: "+str(nearWallCellThicknessRounded))
print("Number of boundary layer points normal to wall: "+str(blThicknessNumPoints))

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
teUpperSpline = gmsh.model.geo.addSpline(teUpperPoints)
leSpline = gmsh.model.geo.addSpline(lePoints)
teLowerSpline = gmsh.model.geo.addSpline(teLowerPoints)
gmsh.model.geo.mesh.setTransfiniteCurve(teUpperSpline, blTENumPoints, "Progression", blTEProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(leSpline, blLENumPoints, "Bump", blLEBump)
gmsh.model.geo.mesh.setTransfiniteCurve(teLowerSpline, blTENumPoints, "Progression", -blTEProgression)

################################################################################
#   Boundary Layer Geometry
################################################################################

#Create boundary layer points
blTopRightPoint = gmsh.model.geo.addPoint(1, blThickness, 0)
blBottomRightPoint = gmsh.model.geo.addPoint(1, -blThickness, 0)
blTopMidPoint = gmsh.model.geo.addPoint(coords[thickestIndex][0], blThickness, 0)
blBottomMidPoint = gmsh.model.geo.addPoint(coords[thickestIndex][0], -blThickness, 0)
blNosePoints = []
for i in range(thickestIndex, len(coords)-thickestIndex):
    blNosePoints.append(gmsh.model.geo.addPoint(coords[i][0], coords[i][1], 0))

#Create the boundary layer lines
blTopTELine = gmsh.model.geo.addLine(blTopRightPoint, blTopMidPoint)
blBottomTELine = gmsh.model.geo.addLine(blBottomMidPoint, blBottomRightPoint)
blTopRightLine = gmsh.model.geo.addLine(aerofoilPoints[0], blTopRightPoint)
blBottomRightLine = gmsh.model.geo.addLine(blBottomRightPoint, aerofoilPoints[-1])
blTopMidLine = gmsh.model.geo.addLine(blTopMidPoint, aerofoilPoints[thickestIndex])
blBottomMidLine = gmsh.model.geo.addLine(aerofoilPoints[len(coords)-thickestIndex-1], blBottomMidPoint)
blLESpline = gmsh.model.geo.addSpline(blNosePoints)
gmsh.model.geo.dilate([(1, blLESpline)], coords[thickestIndex][0], 0, 0, 1+blThickness/coords[thickestIndex][0], blThickness/coords[thickestIndex][1], 1)
gmsh.model.geo.mesh.setTransfiniteCurve(blTopTELine, blTENumPoints, "Progression", blTEProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(blBottomTELine, blTENumPoints, "Progression", -blTEProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(blTopRightLine, blThicknessNumPoints, "Progression", blThicknessProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(blBottomRightLine, blThicknessNumPoints, "Progression", -blThicknessProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(blTopMidLine, blThicknessNumPoints, "Progression", -blThicknessProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(blBottomMidLine, blThicknessNumPoints, "Progression", blThicknessProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(blLESpline, blLENumPoints, "Bump", blLEBump)

#Create boundary layer curve loops
blTopTELoop = gmsh.model.geo.addCurveLoop([blTopTELine, blTopMidLine, -teUpperSpline, blTopRightLine])
blLELoop = gmsh.model.geo.addCurveLoop([blLESpline, -blBottomMidLine, -leSpline, -blTopMidLine])
blBottomTELoop = gmsh.model.geo.addCurveLoop([-teLowerSpline, blBottomMidLine, blBottomTELine, blBottomRightLine])

#Create boundary layer surfaces
blTopTESurface = gmsh.model.geo.addPlaneSurface([blTopTELoop])
blLESurface = gmsh.model.geo.addPlaneSurface([blLELoop])
blBottomTESurface = gmsh.model.geo.addPlaneSurface([blBottomTELoop])
gmsh.model.geo.mesh.setTransfiniteSurface(blTopTESurface)
gmsh.model.geo.mesh.setTransfiniteSurface(blLESurface)
gmsh.model.geo.mesh.setTransfiniteSurface(blBottomTESurface)
gmsh.model.geo.mesh.setRecombine(2, blTopTESurface)
gmsh.model.geo.mesh.setRecombine(2, blLESurface)
gmsh.model.geo.mesh.setRecombine(2, blBottomTESurface)

################################################################################
#   Wake Geometry
################################################################################

#Create points for wake
wakeTopLeftPoint = blTopRightPoint
wakeBottomLeftPoint = blBottomRightPoint
wakeTopRightPoint = gmsh.model.geo.addPoint(wakeLength, wakeEndThickness, 0)
wakeBottomRightPoint = gmsh.model.geo.addPoint(wakeLength, -wakeEndThickness, 0)

#Create lines for wake
wakeTopLine = gmsh.model.geo.addLine(wakeTopRightPoint, wakeTopLeftPoint)
wakeBottomLine = gmsh.model.geo.addLine(wakeBottomLeftPoint, wakeBottomRightPoint)
wakeRightLine = gmsh.model.geo.addLine(wakeBottomRightPoint, wakeTopRightPoint)
gmsh.model.geo.mesh.setTransfiniteCurve(wakeTopLine, wakeNumPoints, "Progression", -wakeProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(wakeBottomLine, wakeNumPoints, "Progression", wakeProgression)
gmsh.model.geo.mesh.setTransfiniteCurve(wakeRightLine, blThicknessNumPoints*2-1)

#Create curveloop for wake
wakeLoop = gmsh.model.geo.addCurveLoop([wakeTopLine, -blTopRightLine, -blBottomRightLine, wakeBottomLine, wakeRightLine])

#Create surface for wake
wakeSurface = gmsh.model.geo.addPlaneSurface([wakeLoop])
gmsh.model.geo.mesh.setTransfiniteSurface(wakeSurface, "Left", [wakeTopRightPoint, wakeTopLeftPoint, wakeBottomLeftPoint, wakeBottomRightPoint])
gmsh.model.geo.mesh.setRecombine(2, wakeSurface)

################################################################################
#   Far Field Geometry
################################################################################

#Create points for the far field
ffTopLeftPoint = gmsh.model.geo.addPoint(-ffHorzOffset, ffHeight/2, 0)
ffTopRightPoint = gmsh.model.geo.addPoint(wakeLength, ffHeight/2, 0)
ffBottomLeftPoint = gmsh.model.geo.addPoint(-ffHorzOffset, -ffHeight/2, 0)
ffBottomRightPoint = gmsh.model.geo.addPoint(wakeLength, -ffHeight/2, 0)

#Create lines for the far field
ffTopLine = gmsh.model.geo.addLine(ffTopRightPoint, ffTopLeftPoint)
ffLeftLine = gmsh.model.geo.addLine(ffTopLeftPoint, ffBottomLeftPoint)
ffBottomLine = gmsh.model.geo.addLine(ffBottomLeftPoint, ffBottomRightPoint)
ffTopRightLine = gmsh.model.geo.addLine(wakeTopRightPoint, ffTopRightPoint)
ffBottomRightLine = gmsh.model.geo.addLine(ffBottomRightPoint, wakeBottomRightPoint)
gmsh.model.geo.mesh.setTransfiniteCurve(ffTopLine, round(ffPointsPerLength*(ffHorzOffset+wakeLength)))
gmsh.model.geo.mesh.setTransfiniteCurve(ffLeftLine, round(ffPointsPerLength*ffHeight))
gmsh.model.geo.mesh.setTransfiniteCurve(ffBottomLine, round(ffPointsPerLength*(ffHorzOffset+wakeLength)))
gmsh.model.geo.mesh.setTransfiniteCurve(ffBottomRightLine, round(ffPointsPerLength*(ffHeight/2-wakeEndThickness/2)))
gmsh.model.geo.mesh.setTransfiniteCurve(ffTopRightLine, round(ffPointsPerLength*(ffHeight/2-wakeEndThickness/2)))

#Create the curve loop for the far field
farFieldLoop = gmsh.model.geo.addCurveLoop([ffTopLine, ffLeftLine, ffBottomLine, ffBottomRightLine, -wakeBottomLine, -blBottomTELine, -blLESpline, -blTopTELine, -wakeTopLine, ffTopRightLine])

#Create the surface
farFieldSurface = gmsh.model.geo.addPlaneSurface([farFieldLoop])

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
farFieldGroup = gmsh.model.addPhysicalGroup(1, [ffTopLine, ffLeftLine, ffBottomLine, ffBottomRightLine, wakeRightLine, ffTopRightLine])
gmsh.model.setPhysicalName(1, farFieldGroup, "farfieldboundary")
fluidGroup = gmsh.model.addPhysicalGroup(2, [farFieldSurface, blLESurface, blTopTESurface, blBottomTESurface, wakeSurface])
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
gmsh.option.setNumber("Mesh.ElementOrder", 3)
gmsh.option.setNumber("Mesh.HighOrderOptimize", 3)
gmsh.option.setNumber("Mesh.NumSubEdges", 10)
gmsh.model.mesh.generate(2)
gmsh.write(aerofoilName+".msh")

#Visualize model
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

#Finalize
gmsh.finalize()