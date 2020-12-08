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
aerofoilGroup = gmsh.model.addPhysicalGroup(2, [aerofoilLoop])
gmsh.model.setPhysicalName(2, aerofoilGroup, "Aerofoil Surface")

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
volumeGroup = gmsh.model.addPhysicalGroup(2, [volumeLoop])
gmsh.model.setPhysicalName(2, volumeGroup, "Volume Boundary")

#Create the surface
surface = gmsh.model.geo.addPlaneSurface([aerofoilLoop, volumeLoop])

#Synchronize CAD entities with the Gmsh model
gmsh.model.geo.synchronize()

#Generate and save mesh
gmsh.model.mesh.generate(2)
gmsh.write("NACA0012.msh")

#Visualize model
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

#Finalize
gmsh.finalize()