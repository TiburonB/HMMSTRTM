# Atom.py 
# TLB 08/23/2020
# This file contains the Atom object and associated functions
# An atom is supplied all data given on a pdb "ATOM" line
# depending on the purposes of the project, values for each data point can be masked
import numpy as np
from Points import *


# TLB 11/10/2020
# adding superpositioning ability to atomset using Kabsch
# This program utilizes the Kabsch algorthm (code can be found here:)
# https://github.com/charnley/rmsd
import rmsd

class Atom:    
    def __init__ (self, Point, name, occupancy, temp_factor, element, res_num = -1, res_type = "", chain_id = "", atom_num = -1, line = ""):
        self.line = line
        self.point = Point
        self.name = name
        self.occupancy = occupancy
        self.temp_factor = temp_factor
        self.element = element
        # Traditionally mask these (loading from *_maskresprot)
        self.res_type = res_type
        self.chain_id = chain_id
        self.atom_num = atom_num
        self.res_num = res_num

    @classmethod
    def from_line_maskresprot(cls, line, prot_masking = False):
        if line[:5] == "ATOM ":
            #line_type = line[:7].strip()
            #atom_num = int(line[7:12])
            name = line[12:17].strip()
            if prot_masking: # special proton masking
                keep_names = [ "HA", "HA2", "HA3", "HN", "HB", "H", "H1", "H2", "H3" ]
                if name not in keep_names:
                    # names of H's on methyl group get renamed to "HM"
                    met_names = ['H1', 'H2', 'H3', 'HE1', 'HE2', 'HE3', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23', 'HD11', 'HD12', 'HD13', 'HD1', 'HD2', 'HD3', 'HZ1', 'HZ2', 'HZ3', 'HD21', 'HD22', 'HD23']
                    if name in met_names:
                        name = "HM"
                    else:
                        name = "HR"

            #res_type = line[17:21].strip()
            #chain_id = line[21:22].strip()
            #res_num = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            occupancy = float(line[54:60])
            temp_factor = float(line[60:66])
            element = line[66:len(line)].strip()
            point = Point(x,y,z)
            return cls(point, name, occupancy, temp_factor, element)
        else:
            print("Invalid line passed into atom constructor method.")
            return -1

    @classmethod
    def from_line(cls, line):
        if line[:5] == "ATOM ":
            atom_num = int(line[7:12])
            name = line[12:17].strip()
            res_type = line[17:21].strip()
            chain_id = line[21:22].strip()
            res_num = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            occupancy = float(line[54:60])
            temp_factor = float(line[60:66])
            element = line[66:len(line)].strip()
            point = Point(x,y,z)
            #print(chain_id)
            return cls(point, name, occupancy, temp_factor, element, res_num, res_type, chain_id, atom_num, line = line)
        else:
            print("Invalid line passed into atom constructor method.")
            return -1


    def __str__(self):
        s = "ATOM " 
        s += str(self.atom_num).rjust(6)
        s += "  " 
        s += str(self.name).ljust(3)
        s += str(self.res_type).rjust(4)
        s += str(self.chain_id).rjust(2)
        s += str(self.res_num).rjust(4)
        s += str(self.point.x).rjust(11)
        s += str(self.point.y).rjust(8)
        s += str(self.point.z).rjust(8)
        s += str(self.occupancy).rjust(5)
        s += "{:.2f}".format(self.temp_factor).rjust(8)
        s += str(self.element).rjust(12)
        s += "\n"
        return s

    def get_point(self):
        return self.point

    def distance(self, Atom2):
        return self.point.distance(Atom2.point)


HYDROGEN_NAMES = ['H', 'HA', 'HB1', 'HB2', 'HB3', 'HG2', 'HG3', 'HD2', 'HD3', 'HE', 'HH11', 'HH12', 'HH21', 'HH22', 'HD21', 'HD22', 'HG', 'HE21', 'HE22', 'HA2', 'HA3', 'HD1', 'HE1', 'HB', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23', 'HD11', 'HD12', 'HD13', 'HD23', 'HE2', 'HE3', 'HZ1', 'HZ2', 'HZ3', 'HZ', 'HG1', 'HH2', 'HH', 'HG11']

# An AtomSet could be a rotamer, a residue, a whole protein
# naming AtomSet optional
class AtomSet:
    def __init__ (self, atoms, name = ""):
        self.name = name
        self.atoms = atoms

        self.is_restricted = False
        self.n = len(atoms)
    
        self.coords = np.zeros((self.n,3))

        for atom_i in range(0,self.n):
            point = self.atoms[atom_i].point
            self.coords[atom_i][:] = [point.x, point.y, point.z]

    def translate(self, tv):
        self.coords = self.coords + tv
        for atom_i in range(0, self.n):
            point = Point.from_array(self.coords[:][atom_i])
            self.atoms[atom_i].point = point

    def rotate(self, rm):
        for atom_i in range(0, self.n):
            self.coords[atom_i][:] = np.matmul(self.coords[atom_i][:], rm)
            self.atoms[atom_i].point = Point.from_array(self.coords[atom_i][:])

    def superpose(self, cloud):

        cloud_centroid = rmsd.centroid(cloud.coords)
        master_centroid = rmsd.centroid(self.coords)

        # center cloud and master at [0,0,0]
        cloud.translate(-cloud_centroid)
        self.translate(-master_centroid)


        # rotate cloud to fit master 
        rm = np.array(rmsd.kabsch(cloud.coords,self.coords))
        cloud.rotate(rm)

        # center cloud at the center of master
        cloud.translate(master_centroid)
        self.translate(master_centroid)

    @classmethod 
    def from_file (cls, filename, name = ""):
        a = []
        f = open( filename, "r")
        line_num = 0
        for line in f:
            if line[:5] == "ATOM ":
                ta = Atom.from_line(line)
                a.append(ta)
        return cls(a)

    def __str__ (self):
        s = ""
        s += self.name + "\n"
        for atom in self.atoms:
            s += str(atom)
        return s

    def write(self, out_f, chain):
        s = ""
        for atom in self.atoms:
            if atom.chain_id == chain:
                s += str(atom)
        f = open(out_f, 'w')
        f.write(s)
        f.close()
        return


    def restrict_to_atoms(self, Anames):
        self.is_restricted = True
        self.restricted_names = Anames
        new_atoms = self.atoms
        new_coords = self.coords

        i = 0
        while i < len(new_atoms):
            if new_atoms[i].name not in Anames:
                new_atoms.pop(i)
                new_coords = np.delete(new_coords, i, axis = 0)
#                print(len(new_atoms))
                i = i - 1
            i = i + 1

        self.n = len(new_atoms)

        #for atom in new_atoms:
            #print(atom.name)

        #print("new n = " + str(self.n))
        #print("new len(new_atoms) = " + str(len(new_atoms)))
        #print("new np.shape(new_coords) = " + str(np.shape(new_coords)))

        self.atoms = new_atoms
        self.coords = new_coords


    # When we want to restrict an AtomSet to just hydrogens...
    # used in cap_ippds.py
    def restrict_to_hydrogens( self ):
        self.restrict_to_atoms(HYDROGEN_NAMES)

    # create a new object variable (self.iads)
    # generate a np array that is a pairwise distance
    # matrix of all atoms in this AtomSet
    def get_iads(self): 
        print("getting iads for " + self.name + ".\n")
        self.iads = np.zeros((self.n, self.n))
        for i in range(0, self.n):
            for j in range(i+1, self.n):
                self.iads[i][j] = self.atoms[i].distance(self.atoms[j])

    def get_iads_1d(self, chain = None): # get iads in 1 dimension
        self.iads = []
        n2 = 0
        for i in range(0, self.n):
            for j in range(i+1, self.n):
                if chain is not None and self.atoms[i].chain_id == chain and self.atoms[j].chain_id == chain:
                    self.iads.append(self.atoms[i].distance(self.atoms[j]))
                    n2 += 1

        print(str(math.sqrt(n2)))
        return self.iads

    # function returns new atomset that is a subset of the current set
    # indices for the subsection can be defined, name changed accordingly:
    # new_name = <old_name>_subsect_<lower_bound>:<upper_bound>
    # where lower_bound and upper_bound mark the subset of the original 
    def subset_atoms(self, l = 0, u = -1):
        if u == -1:
            u = self.n
        new_name = self.name + "_subsect_" + str(l) + ":" + str(u)
        new_set = AtomSet(self.atoms[l:u], new_name)
        new_set.is_restricted = self.is_restricted

        if self.is_restricted:
            new_set.restricted_names = self.restricted_names 

        return new_set










