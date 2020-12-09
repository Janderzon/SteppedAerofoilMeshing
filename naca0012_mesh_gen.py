import naca_4_series_points
import gmsh
import sys

#Parameters
n = 50     #Number of panels for upper and lower surface
lc = 1e-1   #Target mesh size close to points
s = "sin"   #Spacing type for aerofoil points

#Get the coordinates of the aerofoil
coords = naca_4_series_points.points(0, 0, 12, n, s)

#Initialize Gmsh and name the model
gmsh.initialize()
gmsh.model.add("NACA 0012")

#Create the aerofoil points
upperPoints = []
lowerPoints = []
for i in range(0, n+1):
    upperPoints.append(gmsh.model.geo.addPoint(coords[0][i], coords[1][i], 0, lc))
    if i>0:
        lowerPoints.append(gmsh.model.geo.addPoint(coords[2][i], coords[3][i], 0, lc))

#Create the aerofoil lines
upperLines = []
lowerLines = []
for i in range(0, n):
    upperLines.append(gmsh.model.geo.addLine(upperPoints[i], upperPoints[i+1]))
    if i>0:
        lowerLines.append(-gmsh.model.geo.addLine(lowerPoints[i-1], lowerPoints[i]))
    else:
        lowerLines.append(-gmsh.model.geo.addLine(upperPoints[i], lowerPoints[i]))
upperLines.append(gmsh.model.geo.addLine(upperPoints[n], lowerPoints[n-1]))
lowerLines.reverse()

#Create the aerofoil curve loop
aerofoilLoop = gmsh.model.geo.addCurveLoop(lowerLines+upperLines)

#Create boundary layer points
upperBLPoints = []
lowerBLPoints = []
for i in range(0, n+1):
    upperBLPoints.append((0,gmsh.model.geo.addPoint(coords[0][i], coords[1][i], 0, lc)))
    if i>0:
        lowerBLPoints.append((0,gmsh.model.geo.addPoint(coords[2][i], coords[3][i], 0, lc)))

#Scale boundary layer points
gmsh.model.geo.dilate(upperBLPoints+lowerBLPoints, 0.25, 0, 0, 1.5, 2, 1)

#Create the boundary layer lines
upperBLLines = []
lowerBLLines = []
for i in range(0, n):
    upperBLLines.append(gmsh.model.geo.addLine(upperBLPoints[i][1], upperBLPoints[i+1][1]))
    if i>0:
        lowerBLLines.append(-gmsh.model.geo.addLine(lowerBLPoints[i-1][1], lowerBLPoints[i][1]))
    else:
        lowerBLLines.append(-gmsh.model.geo.addLine(upperBLPoints[i][1], lowerBLPoints[i][1]))
upperBLLines.append(gmsh.model.geo.addLine(upperBLPoints[n][1], lowerBLPoints[n-1][1]))
lowerBLLines.reverse()

#Create the boundary layer curve loop
boundaryLayerLoop = gmsh.model.geo.addCurveLoop(lowerBLLines+upperBLLines)

#Create lines to divide up boundary layer into quadrilaterals
upperBLDivideLines = []
lowerBLDivideLines = []
for i in range(0, n+1):
    upperBLDivideLines.append(gmsh.model.geo.addLine(upperBLPoints[i][1], upperPoints[i]))
    if i<n:
        lowerBLDivideLines.append(gmsh.model.geo.addLine(lowerBLPoints[i][1], lowerPoints[i]))
    #else:
        #lowerBLDivideLines.append(-gmsh.model.geo.addLine(upperBLPoints[i][1], lowerBLPoints[i][1]))
#upperBLDivideLines.append(gmsh.model.geo.addLine(upperBLPoints[n][1], lowerBLPoints[n-1][1]))
#lowerBLDivideLines.reverse()

#Create points for the volume boundary
topLeft = gmsh.model.geo.addPoint(-1, 1, 0, lc)
topRight = gmsh.model.geo.addPoint(2, 1, 0, lc)
bottomLeft = gmsh.model.geo.addPoint(-1, -1, 0, lc)
bottomRight = gmsh.model.geo.addPoint(2, -1, 0, lc)

#Create lines for the volume boundary
top = gmsh.model.geo.addLine(topLeft, topRight)
right = gmsh.model.geo.addLine(topRight, bottomRight)
bottom = gmsh.model.geo.addLine(bottomRight, bottomLeft)
left = gmsh.model.geo.addLine(bottomLeft, topLeft)

#Create the curve loop for the volume boundary
volumeLoop = gmsh.model.geo.addCurveLoop([top, right, bottom, left])

#Create the surface
surface = gmsh.model.geo.addPlaneSurface([boundaryLayerLoop, volumeLoop])

#Synchronize CAD entities with the Gmsh model
gmsh.model.geo.synchronize()

#Create physical groups
aerofoilGroup = gmsh.model.addPhysicalGroup(1, [aerofoilLoop])
gmsh.model.setPhysicalName(1, aerofoilGroup, "Aerofoil Surface")
volumeGroup = gmsh.model.addPhysicalGroup(1, [volumeLoop])
gmsh.model.setPhysicalName(1, volumeGroup, "Volume Boundary")

#Generate and save mesh
gmsh.model.mesh.generate(2)
gmsh.write("NACA0012.msh")

#Visualize model
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

#Finalize
gmsh.finalize()