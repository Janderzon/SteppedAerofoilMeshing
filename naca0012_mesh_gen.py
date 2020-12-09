import naca_4_series_points
import gmsh
import sys
import aerofoil_scale

#Parameters
n = 25                          #Number of panels for each of the upper and lower surfaces of the aerofoil
spacing = "sin"                 #Spacing type for aerofoil points
farFieldSize = 10               #Width and height of the far field volume boundary
farFieldMeshSize = 1            #Target mesh size at the far field volume boundary
boundaryLayerMeshSize = 1e-1    #Target mesh size at the edge of the boundary layer
boundaryLayerThickness = 0.02   #Thickness of the boundary layer
boundaryLayerNumCells = 10      #Number of cells in the boundary layer in the direction normal to the aerofoil surface
boundaryLayerProgression = 1.5  #Growth rate of the cells in the boundary layer

#Get the coordinates of the aerofoil
coords = naca_4_series_points.points(0, 0, 12, n, spacing)

#Initialize Gmsh and name the model
gmsh.initialize()
gmsh.model.add("NACA 0012")

#Create the aerofoil points
aerofoilPoints = []
for i in range(0, len(coords)):
    aerofoilPoints.append(gmsh.model.geo.addPoint(coords[i][0], coords[i][1], 0))

#Create the aerofoil lines
aerofoilLines = []
for i in range(0, len(aerofoilPoints)):
    aerofoilLines.append(gmsh.model.geo.addLine(aerofoilPoints[i-1], aerofoilPoints[i]))
    gmsh.model.geo.mesh.setTransfiniteCurve(aerofoilLines[i],2)

#Create the aerofoil curve loop
aerofoilLoop = gmsh.model.geo.addCurveLoop(aerofoilLines)

#Create boundary layer points
boundaryLayerCoords = aerofoil_scale.boundaryLayerScale(coords,boundaryLayerThickness)
boundaryLayerPoints = []
for i in range(0, len(boundaryLayerCoords)):
    boundaryLayerPoints.append(gmsh.model.geo.addPoint(boundaryLayerCoords[i][0], boundaryLayerCoords[i][1], 0, boundaryLayerMeshSize))

#Create the boundary layer lines
boundaryLayerLines = []
for i in range(0, len(boundaryLayerPoints)):
    boundaryLayerLines.append(gmsh.model.geo.addLine(boundaryLayerPoints[i-1], boundaryLayerPoints[i]))
    gmsh.model.geo.mesh.setTransfiniteCurve(boundaryLayerLines[i],2)

#Create the boundary layer curve loop
boundaryLayerLoop = gmsh.model.geo.addCurveLoop(boundaryLayerLines)

#Create lines to divide up boundary layer into quadrilaterals
boundaryLayerDividerLines = []
for i in range(0, len(boundaryLayerPoints)):
    boundaryLayerDividerLines.append(gmsh.model.geo.addLine(aerofoilPoints[i], boundaryLayerPoints[i]))
    gmsh.model.geo.mesh.setTransfiniteCurve(boundaryLayerDividerLines[i],boundaryLayerNumCells,"Progression",boundaryLayerProgression)

#Create boundary layer quadrilateral curve loops
boundaryLayerQuadLoops = []
for i in range(0, len(boundaryLayerPoints)):
    boundaryLayerQuadLoops.append(gmsh.model.geo.addCurveLoop([boundaryLayerLines[i],-boundaryLayerDividerLines[i],-aerofoilLines[i],boundaryLayerDividerLines[i-1]]))

#Create boundary layer quadrilateral surfaces
boundaryLayerQuadSurfaces = []
for i in range(0, len(boundaryLayerPoints)):
    boundaryLayerQuadSurfaces.append(gmsh.model.geo.addPlaneSurface([boundaryLayerQuadLoops[i]]))
    gmsh.model.geo.mesh.setTransfiniteSurface(boundaryLayerQuadSurfaces[i])
    gmsh.model.geo.mesh.setRecombine(2,boundaryLayerQuadSurfaces[i])

#Create points for the volume boundary
topLeft = gmsh.model.geo.addPoint(-farFieldSize/2+0.5, farFieldSize/2, 0, farFieldMeshSize)
topRight = gmsh.model.geo.addPoint(farFieldSize/2+0.5, farFieldSize/2, 0, farFieldMeshSize)
bottomLeft = gmsh.model.geo.addPoint(-farFieldSize/2+0.5, -farFieldSize/2, 0, farFieldMeshSize)
bottomRight = gmsh.model.geo.addPoint(farFieldSize/2+0.5, -farFieldSize/2, 0, farFieldMeshSize)

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