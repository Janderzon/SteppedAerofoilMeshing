import math

def boundaryLayerScale(coords, size):
    newCoords = []
    for i in range(0, len(coords)):
        #Create vectors
        vector1 = []
        vector2 = []
        if i<len(coords)-1:
            vector1 = [coords[i][0]-coords[i-1][0], coords[i][1]-coords[i-1][1]]
            vector2 = [coords[i+1][0]-coords[i][0], coords[i+1][1]-coords[i][1]]
        else:
            vector1 = [coords[i][0]-coords[i-1][0], coords[i][1]-coords[i-1][1]]
            vector2 = [coords[0][0]-coords[i][0], coords[0][1]-coords[i][1]]

        #Normalise vectors
        vector1Mag = math.sqrt(vector1[0]**2+vector1[1]**2)
        vector2Mag = math.sqrt(vector2[0]**2+vector2[1]**2)
        for j in range(0, 2):
            vector1[j] /= vector1Mag
            vector2[j] /= vector2Mag
        
        #Sum vectors
        totalVector = [vector1[0]+vector2[0], vector1[1]+vector2[1]]

        #Find the normal to the surface
        normal = [totalVector[1], -totalVector[0]]

        #Normalise the normal
        normalMag = math.sqrt(normal[0]**2+normal[1]**2)
        normal[0] /= normalMag
        normal[1] /= normalMag

        #Translate the coordinates along the normal
        newCoords.append((coords[i][0]+normal[0]*size, coords[i][1]+normal[1]*size))
        
    return newCoords