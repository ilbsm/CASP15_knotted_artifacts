import numpy as np
from ast import literal_eval
from scipy.spatial import distance
from collections import defaultdict, Counter
from icecream import ic
from jones import jones
from tqdm import tqdm
#from topoly import Graph

def projection2gcode(projection):
    stacked = np.outer(projection, np.roll(projection, -1, axis=1))
    ic(stacked.shape)
    min_x_seg = np.min.outer(projection[0], np.roll(projection[0], -1, axis=0))
    max_x_seg = np.max.outer(projection[0], np.roll(projection[0], -1, axis=0))
    ic(min_seg)

    segments = (np.roll(projection, -1, axis=1) - projection)[:,:-1]
    dx,dy,dz = segments
    Sx = (np.roll(projection[0], -1, axis=0) + projection[0])[:-1]
    Sy = (np.roll(projection[1], -1, axis=0) + projection[1])[:-1]
    # y=ax+b -- equation of line defined by a segment
    a = dy/dx
    b = (Sy-a*Sx)/2
    da = np.subtract.outer(a, a)
    db = np.subtract.outer(b, b)
    inv_a = 1/a
    b_inv_a = b*inv_a
    d_inv_a = np.subtract.outer(inv_a, inv_a)
    db_inv_a = np.subtract.outer(b_inv_a, b_inv_a)
    x_intersection = -db/da
    y_intersection = db_inv_a/d_inv_a
    seg_div_x = (x_intersection-projection[0,:-1])/dx
    seg_div_y = (y_intersection-projection[1,:-1])/dy
    where_cross_x = np.logical_and(0<seg_div_x, seg_div_x<1)
    where_cross_y = np.logical_and(0<seg_div_y, seg_div_y<1)
    where_cross = np.logical_and(where_cross_x, where_cross_y)
    where_cross = np.argwhere(np.logical_and(where_cross, where_cross.T))
    where_cross = [(ndx1,ndx2) for ndx1,ndx2 in where_cross if ndx1+1<ndx2]
    crossings = []
    for ndx1, ndx2 in where_cross:
        is_overcrossing = get_is_overcrossing(projection, ndx1, ndx2, seg_div_x)
        sign = get_crossing_sign(segments, ndx1, ndx2, is_overcrossing)
        s1 = round(seg_div_x[ndx1, ndx2],3)
        s2 = round(seg_div_x[ndx2, ndx1],3)
        crossings.append((ndx1+s1, ndx2+s2, -1+2*is_overcrossing, -1+2*sign))
    gcode = crossings2gcode(sorted(crossings))
    return gcode

def get_crossing_sign(segments, ndx1, ndx2, is_overcrossing):
    cross_product = np.cross(segments[:2,ndx1], segments[:2,ndx2], axisa=0, axisb=0)
    if cross_product == 0:
        raise ValueError
    sign = is_overcrossing == (cross_product>0)
    return sign

def get_is_overcrossing(projection, ndx1, ndx2, seg_div):
    z1 = projection[2,ndx1]
    z2 = projection[2,ndx1+1]
    z3 = projection[2,ndx2]
    z4 = projection[2,ndx2+1]
    z_cross_12 = seg_div[ndx1,ndx2] * (z2-z1) + z1
    z_cross_34 = seg_div[ndx1,ndx2] * (z4-z3) + z3
    is_overcrossing = z_cross_34 < z_cross_12
    return is_overcrossing

def crossings2gcode(crossings):
    symbols = []
    indices = []
    is_over = []
    sign_dict = {}
    i=1
    for i, (ndx1, ndx2, io, sign) in enumerate(crossings, start=1):
        symbols += [i,i]
        indices += [ndx1, ndx2]
        is_over += [io, -io]
        sign_dict[i] = sign
    order = np.argsort(indices)
    symbols = [0] + list(np.array(symbols)[order]) + [i+1]
    is_over = [0] + list(np.array(is_over)[order]) + [0]
    if is_over[1] == -1:
        is_over = [-x for x in is_over]
    sign_dict[0] = 0
    sign_dict[i+1] = 0
    gcode = ([symbols], [is_over], sign_dict)
    return gcode

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
    coords = np.einsum('ijn,jk->ikn',BA,coords)
    return coords

def kmt(coords):
    a = coords.T.tolist()
    g = Graph(a)
    g.reduce()
    new_coords = []
    for k in sorted(g.coordinates.keys()):
        new_coords.append(g.coordinates[k])
    return np.array(new_coords).T

def group_gcodes(gcodes):
    group = Counter()
    for gcode in gcodes:
        gcode = str(gcode).replace(' ','')
        group[gcode] += 1
    all_counts = 0
    res = []
    for _, count in group.most_common():
        all_counts += count
    for gcode, count in group.most_common():
        res.append((count/all_counts, literal_eval(gcode)))
    return res

def gcodes2jones(gcode_groups):
    all_counts = 0
    results = Counter()
    for counts, data in gcode_groups:
        all_counts += counts
        gcode, isover, sign_dict = data
        res = jones(gcode, isover, sign_dict, is_knotoid=True, use_canonical=False, 
                    show_results=False, debug=False, to_write=False)
        results[str(res).replace(' ','')] += counts
    results = [(c/all_counts,*literal_eval(k)) for k,c in results.most_common()]
    return results # list of (weight, polynomial, lowest_exponent)
    
def get_mean_jones(polynomials): # output of gcodes2jones as input
    lowest_exp = np.inf
    highest_exp = -np.inf
    for _, poly, min_exp in polynomials:
        lowest_exp = min(lowest_exp, min_exp)
        highest_exp = max(min_exp + (len(poly)-1)/2, highest_exp)
    final_poly = [0]*int((highest_exp-lowest_exp)*2+1)
    for weight, poly, min_exp in polynomials:
        first_index = (min_exp - lowest_exp)*2
        assert first_index.is_integer()
        first_index = int(first_index)
        for i, p in enumerate(poly):
            final_poly[i+first_index] += p*weight
    return final_poly, lowest_exp

def compare_target_model(polynomials1, polynomials2): # outputs of gcodes2jones as input
    polys1_dict = defaultdict(int)
    polys2_dict = defaultdict(int)
    keys = set([])
    for weight, poly, lowest_exp in polynomials1:
        key = str((poly,lowest_exp)).replace(' ','')
        polys1_dict[key] = weight
        keys.add(key)
    for weight, poly, lowest_exp in polynomials2:
        key = str((poly,lowest_exp)).replace(' ','')
        polys2_dict[key] = weight
        keys.add(key)
    similarity = 0
    for key in keys:
        similarity += min(polys1_dict[key], polys2_dict[key])
    return similarity
   
def coords2jones(coords, num_rotations=1):
    projections = rotate_randomly(coords, num_rotations)
    projections = coords.reshape(*coords.shape,1)
    gcodes = []
    for i in tqdm(range(projections.shape[-1])):
        projection = projections[:,:,i]
#        projection = kmt(projection)
        gcode = projection2gcode(projection)
        gcodes.append(gcode)
    gcode_groups = group_gcodes(gcodes)
    for gcode in gcode_groups:
        print(gcode)
#    print(len(gcode_groups))
    jones_list = gcodes2jones(gcode_groups)
    return jones_list
    

if __name__ == '__main__':
    #coords, ndxs = load_nxyz('R1117TS416_3.xyz')
    #coords, ndxs = load_xyz('31.xyz')
    coords, ndxs = load_xyz('31.xyz')
    jones_list1 = coords2jones(coords)
    #jones_list2 = get_jones('89.xyz')
    print(jones_list1)
    #print(jones_list2)
    print(get_mean_jones(jones_list1))
    #print(get_mean_jones(jones_list2))
    #print(compare_target_model(jones_list1, jones_list2))
    




