import numpy as np
from icecream import ic

def get_diagrams(projections):
    delta = np.zeros((projections.shape[0]-1, projections.shape[1]-1, *projections.shape[2:]))
    delta[0] = (np.roll(projections[0], -1, axis=0) - projections[0])[:-1]
    delta[1] = (np.roll(projections[1], -1, axis=0) - projections[1])[:-1]
    dx = delta[0]
    dy = delta[1]
    Sx = (np.roll(projections[0], -1, axis=0) + projections[0])[:-1]
    Sy = (np.roll(projections[1], -1, axis=0) + projections[1])[:-1]
    # y=ax+b -- equation of line defined by a segment
    a = dy/dx
    b = (Sy - a*Sx)/2
    _, num_atoms, num_projections = projections.shape
    x_intersections = np.zeros((num_atoms-1, num_atoms-1, num_projections))
    y_intersections = np.zeros((num_atoms-1, num_atoms-1, num_projections))
    crossings = np.zeros((num_atoms-1, num_atoms-1, num_projections))
    # vectorization?
    for ndx1 in range(num_atoms-1):
        for ndx2 in range(ndx1+2, num_atoms-1):
            x_intersection = -(b[ndx2]-b[ndx1]) / (a[ndx2]-a[ndx1])
            y_intersection = x_intersection * a[ndx1] + b[ndx1]
            # checking if lines' intersection is segments' intersection
            is_crossings = is_intersection_local(projections, ndx1, ndx2,
                                                 x_intersection, y_intersection)
            # checking crossings' signs
            if np.any(is_crossings):
                is_overcrossing = get_is_overcrossing(projections, ndx1, ndx2, x_intersection)
                signs = get_crossing_sign(delta, ndx1, ndx2, is_overcrossing)
                crossings[ndx1,ndx2] = np.where(is_crossings, signs, 0)
    ic(np.argwhere(crossings!=0))
    return crossings
 
def get_crossing_sign(delta, ndx1, ndx2, is_overcrossing):
    cross_product = np.cross(delta[:,ndx1], delta[:,ndx2], axisa=0, axisb=0)
    if np.any(cross_product == 0):
        raise ValueError
    signs = np.where(np.equal(is_overcrossing, cross_product>0), 1, -1)
    return signs

def get_is_overcrossing(projections, ndx1, ndx2, x_intersection):
    x1 = projections[0,ndx1]
    x2 = projections[0,ndx1+1]
    z1 = projections[2,ndx1]
    z2 = projections[2,ndx1+1]
    x3 = projections[0,ndx2]
    x4 = projections[0,ndx2+1]
    z3 = projections[2,ndx2]
    z4 = projections[2,ndx2+1]
    z_cross_12 = (z2*(x_intersection-x1) + z1*(x2-x_intersection))/(x2-x1)
    z_cross_34 = (z4*(x_intersection-x3) + z3*(x4-x_intersection))/(x4-x3)
    is_overcrossing = np.greater(z_cross_34, z_cross_12)
    return is_overcrossing

def is_intersection_local(projections, ndx1, ndx2, x0, y0):
    x1_boundary = np.sort(np.stack((projections[0,ndx1], projections[0,ndx1+1]), axis=1),
                         axis=1)
    y1_boundary = np.sort(np.stack((projections[1,ndx1], projections[1,ndx1+1]), axis=1),
                         axis=1)
    x2_boundary = np.sort(np.stack((projections[0,ndx2], projections[0,ndx2+1]), axis=1),
                         axis=1)
    y2_boundary = np.sort(np.stack((projections[1,ndx2], projections[1,ndx2+1]), axis=1),
                         axis=1)
    # checking if intersection between lines is between endpoints of segments
    bool1 = np.greater(x1_boundary[:,1], x0[:])
    bool2 = np.greater(x0[:], x1_boundary[:,0])
    bool3 = np.greater(x2_boundary[:,1], x0[:])
    bool4 = np.greater(x0[:], x2_boundary[:,0])
    bool5 = np.greater(y1_boundary[:,1], y0[:])
    bool6 = np.greater(y0[:], y1_boundary[:,0])
    bool7 = np.greater(y2_boundary[:,1], y0[:])
    bool8 = np.greater(y0[:], y2_boundary[:,0])
    is_crossings = np.logical_and(np.logical_and(np.equal(bool1, bool2), np.equal(bool3, bool4)),
                                  np.logical_and(np.equal(bool5, bool6), np.equal(bool7, bool8)))
    return is_crossings

def load_nxyz(infile):
    with open(infile, 'r') as f:
        ndxs = []
        coords = []
        for line in f.readlines():
            ndx, x, y, z = line.strip().split()
            ndx = int(ndx)
            ndxs.append(ndx)
            coords.append([float(k) for k in [x,y,z]])
    return np.array(coords).T, ndxs
    # shape(coords) = (3, chain_length)

def load_xyz(infile):
    with open(infile, 'r') as f:
        ndxs = []
        coords = []
        for ndx, line in enumerate(f.readlines()):
            x, y, z = line.strip().split()
            ndxs.append(ndx+1)
            coords.append([float(k) for k in [x,y,z]])
    return np.array(coords).T, ndxs
    # shape(coords) = (3, chain_length)

def rotate_randomly(coords, num_rotations=2):
    # rotated by A(theta)B(phi)
    phi = 2*np.pi*np.random.rand(num_rotations)
    theta = np.pi*np.random.rand(num_rotations)
    # A(-theta)B(-phi)
    # TRANSPOSED, BUT DOESN'T LOOK LIKE (bracket order)
    # transposition handled by np.einsum
    BA = np.array([[ np.cos(phi),         -np.cos(theta)*np.sin(phi),  np.sin(theta)*np.sin(phi) ],
                   [ np.sin(phi),          np.cos(theta)*np.cos(phi), -np.cos(phi)*np.sin(theta) ],
                   [ np.zeros_like(theta), np.sin(theta),              np.cos(theta)             ]])
    return np.einsum('ijn,jk->ikn',BA,coords)

#for debugging
def calc_lengths(coords):
    return np.sqrt(np.sum(np.square(coords - np.roll(coords, -1, axis=1)), axis=0))

if __name__ == '__main__':
    #coords, ndxs = load_nxyz('R1117TS416_3.xyz')
    coords, ndxs = load_xyz('89.xyz')
    # possibly chain reduction (like KMT) here
    num_rotations = 6
    rotated_coords = rotate_randomly(coords, num_rotations)
    ic(rotated_coords.shape)
    get_diagrams(rotated_coords)

