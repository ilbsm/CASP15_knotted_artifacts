import math
import os

class TclFileWriter():
    def __init__(self, last_residue):
        self.name = 'show_ends.tcl'
        self.last = last_residue
        self.tcl = [
            "# VMD for LINUXAMD64, version 1.9.4a57 (April 27, 2022)",
            "# Log file '/home/bagren/lab/CASP15_RNA/show_ends.tcl', created by user bagren",
            "mol modstyle 0 0 Licorice 0.300000 12.000000 12.000000",
            "mol modselect 0 0 backbone and name C2' C3' C4' C5' O3' O4' O5' P",
            "mol modcolor 0 0 Index",
            "mol addrep 0",
            "mol modselect 1 0 resid 1",
            "mol modstyle 1 0 Surf 1.400000 0.000000",
            "mol modcolor 1 0 ColorID 1",
            "mol addrep 0",
            "mol modselect 2 0 resid {:d}".format(self.last),
            "mol modstyle 2 0 Surf 1.400000 0.000000",
            "mol modcolor 2 0 ColorID 0",
            "# VMD for LINUXAMD64, version 1.9.4a57 (April 27, 2022)",
            "# end of log file."]

    def __enter__(self):
        with open(self.name, 'w') as f:
            f.write('\n'.join(self.tcl))

    def __exit__(self, *args):
        os.remove(self.name)

def calc_dist(nxyz1, nxyz2):
    _, x1, y1, z1 = nxyz1.strip().split()
    _, x2, y2, z2 = nxyz2.strip().split()
    x1,x2,y1,y2,z1,z2 = [float(k) for k in [x1,x2,y1,y2,z1,z2]]
    return math.sqrt((x1-x2)**2 +(y1-y2)**2 +(z1-z2)**2)

if __name__ == '__main__':
    files = []
    with open('knots_to_verify.txt', 'r') as f:
        for line in f.readlines():
            infile = line.strip().split('.')[0]
            files.append(infile)

    for infile in files:
        with open('xyz_models/{}.xyz'.format(infile), 'r') as f:
            beg, *_, end = f.readlines()
        with open('../pdb_models/{}.pdb'.format(infile), 'r') as f:
            for line in f.readlines():
                if line[0:4] == 'ATOM':
                    ndx = line.split()[5]
        dist = calc_dist(beg, end)
        print('  file:{}\n  distance between endpoints in angstroms: {:.2f}\n  last residue index: {}'.format(infile, dist, ndx))
        x = input('Show in vmd? (Pass any character to omit this structure)')
        if x:
            continue
        else:
            with TclFileWriter(int(ndx)) as g:
                os.system('vmd ../pdb_models/{}.pdb -e show_ends.tcl'.format(infile))

