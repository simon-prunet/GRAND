import numpy as np

def LinePlaneCollision(planeNormal, planePoint, rayDirection, rayPoint, epsilon=1e-6):
    ndotu = planeNormal.dot(rayDirection)
    if abs(ndotu) < epsilon:
        raise RuntimeError("no intersection or line is within plane")
    w = rayPoint - planePoint
    si = -planeNormal.dot(w) / ndotu
    Psi = w + si * rayDirection + planePoint
    return Psi

def Simplerot(zen_deg):
    zen_rad = np.deg2rad(zen_deg)
    c = np.cos(zen_rad)
    s = np.sin(zen_rad)
    mat = np.array([[c,0.0,-s],[0.0,1.0,0.0],[s,0,c]])
    return mat

