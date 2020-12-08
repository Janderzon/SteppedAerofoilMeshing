import math

def points(m, p, t, n, s=""):
    #Convert inputs to percentages
    m = m*0.01
    p = p*0.1
    t = t*0.01

    #Discretise along the chord
    x = []
    for i in range(0, n+1):
        if s=="-sin":
            x.append(math.sin(math.pi*i/(n*2)))
        elif s=="sin":
            x.append(1-math.cos(math.pi*i/(n*2)))
        elif s=="cos":
            x.append(0.5-0.5*math.cos(math.pi*i/n))
        else:
            x.append(i/n)

    #Get the camber at each point
    yc = []
    for i in x:
        if i<p:
            yc.append((m/(p**2))*(2*p*i-i**2))
        else:
            yc.append((m/((1-p)**2))*((1-2*p)+2*p*i-i**2))

    #Get the thickness above and below the mean line
    yt = []
    for i in x:
        yt.append((t/0.2)*(0.2969*math.sqrt(i)-0.1260*i-0.3516*i**2+0.2843*i**3-0.1015*i**4))

    #Get theta
    theta = []
    for i in x:
        if i<p:
            theta.append(math.atan((m/(p**2))*(2*p-2*i)))
        else:
            theta.append(math.atan((m/((1-p)**2))*(2*p-2*i)))

    #Get the upper surface coordinates
    xu = []
    yu = []
    for i in range(0, n+1):
        xu.append(x[i]-yt[i]*math.sin(theta[i]))
        yu.append(yc[i]+yt[i]*math.cos(theta[i]))

    #Get the lower surface coordinates
    xl = []
    yl = []
    for i in range(0, n+1):
        xl.append(x[i]+yt[i]*math.sin(theta[i]))
        yl.append(yc[i]-yt[i]*math.cos(theta[i]))

    return [xu, yu, xl, yl]