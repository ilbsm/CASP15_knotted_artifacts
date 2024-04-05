import os
from tqdm import tqdm

def pdb2xyz(pdb):
    corr_atoms = ["C5'","C4'","C3'","O3'","P","O5'"]
    coords = []
    i = 1
    for line in pdb.split('\n'):
        if line[:4] == 'ATOM':
            line_split = line.split()
            if line_split[2] in corr_atoms:
                x = line[30:38]
                y = line[38:46]
                z = line[46:54]
                x,y,z = [float(k) for k in (x,y,z)]
                coords.append('{:03d} {:8.3f} {:8.3f} {:8.3f}\n'.format(i,x,y,z))
                i += 1
    return ''.join(coords)

if __name__ == '__main__':
    for folder in ['models', 'targets']:
        if not os.path.exists('xyz_{}'.format(folder)):
            os.mkdir('xyz_{}'.format(folder))
        pdbs = sorted(os.listdir('../pdb_{}'.format(folder)))
        for pdb in tqdm(pdbs):
            if pdb[0] == '.': continue
            name = pdb.split('.')[:-1]
            name = '.'.join(name)
            with open('../pdb_{}/{}'.format(folder, pdb), 'r') as f:
                loaded = f.read()
                xyz = pdb2xyz(loaded)
            with open('xyz_{}/{}.xyz'.format(folder, name), 'w') as g:
                g.write(xyz)
