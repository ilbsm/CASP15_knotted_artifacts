import numpy as np

def get_diagrams(projections):
    # y=ax+b -- equation of line defined by a segment
    # z = ax+by+c
    dx = (projections[0] - np.roll(projections[0], -1, axis=0))[:-1]
    dy = (projections[1] - np.roll(projections[1], -1, axis=0))[:-1]
    dz = np.where(np.greater(projections[2], np.roll(projections[2],1,axis=0)), 1, -1)[:-1]
    print(dz)
    Sx = (projections[0] + np.roll(projections[0], -1, axis=0))[:-1]
    Sy = (projections[1] + np.roll(projections[1], -1, axis=0))[:-1]
    a = dy/dx
    b = (Sy - a*Sx)/2
    _, num_atoms, num_projections = projections.shape
    x_intersections = np.zeros((num_atoms-1, num_atoms-1, num_projections))
    y_intersections = np.zeros((num_atoms-1, num_atoms-1, num_projections))
    crossings = np.zeros((num_atoms-1, num_atoms-1, num_projections))
    for ndx1 in range(num_atoms-1):
        for ndx2 in range(ndx1+1, num_atoms-1):
            x_intersection = -(b[ndx1]-b[ndx2]) / (a[ndx1]-a[ndx2])
            y_intersection = x_intersection * a[ndx1] + b[ndx1]
            x_boundary = np.sort(np.stack((projections[0,ndx1], projections[0,ndx2]), axis=1), axis=1)
            y_boundary = np.sort(np.stack((projections[1,ndx1], projections[1,ndx2]), axis=1), axis=1)
            crossings[ndx1,ndx2] = np.where(np.greater(y_boundary[:,1], y_intersection[:]), 
                                            np.where(np.greater_equal(y_intersection[:], y_boundary[:,0]),
                                                     dz[ndx1], 0), 
                                            0)
    print(crossings)
#            crossings[ndx1,ndx2] = 
#            x_intersections[ndx1,ndx2] = -(b[ndx1]-b[ndx2]) / (a[ndx1]-a[ndx2])
#    y_intersections = a*x_intersections + b
    
def load_xyz(infile):
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
    coords, ndxs = load_xyz('R1117TS416_3.xyz')
    # possibly chain reduction (like KMT) here
    coords = coords[:,:10]
    num_rotations = 1
    rotated_coords = rotate_randomly(coords, num_rotations)
    get_diagrams(rotated_coords)

