import numpy as np
import re
import subprocess
import os
import math

def main():

    if not os.path.exists('results'):
        os.makedirs('results')

    '''
    1. pure bcp film
    '''
    # read in data (by self-avoiding random walk)
    data = LammpsData("data_C20n2400_wo_substate.txt")

    # generate a layer
    layer = Layer(data, film_types=[1, 2])

    # add a substrate
    data_w_substrate = layer.add_substrate(type_A = 3, type_B = 4, ratio=0.5)

    # write the new data
    data_w_substrate.write_data("data_tuned.txt", header="lmp data file")

    '''
    2. bilayer film
    '''
    # read in film(s) (already with their own substrates)
    #data_C = LammpsData("data_C20n4800.txt")

    # generate layers
    #C_top = Layer(data_C, film_types=[1, 2], substrate_types=[3, 4], res_film_types=[1, 2], res_bond_types=[1, 2, 3])
    #C_bottom = Layer(data_C, film_types=[1, 2], substrate_types=[3, 4], res_film_types=[3, 4], res_bond_types=[4, 5, 6])

    # layer1 over layer2; use truediv to place one on top of the other
    #ConC = C_top / C_bottom

    # write the new data
    #ConC.write_data("data_tuned.txt", header="lmp data file of C20n4800 on C20n4800")




class LammpsData:
    '''Lammps data'''
    def __init__(self, fname=None):
        if fname is not None:
            self.fname = os.path.join("ingredients", fname)
            self.read_data()

    def read_data(self):

        flag = 0
        f = open(self.fname, "r")
        while flag == 0:
            line = f.readline()
            match = re.search("atoms", line)
            if match:
                self.num_atoms = int(line.split()[0])
            match = re.search("atom types", line)
            if match:
                self.num_atomtypes = int(line.split()[0])
            match = re.search("bonds", line)
            if match:
                self.num_bonds = int(line.split()[0])
            match = re.search("bond types", line)
            if match:
                self.num_bondtypes = int(line.split()[0])
            match = re.search("xlo xhi", line)
            if match:
                self.xlo = float(line.split()[0])
                self.xhi = float(line.split()[1])
            match = re.search("ylo yhi", line)
            if match:
                self.ylo = float(line.split()[0])
                self.yhi = float(line.split()[1])
            match = re.search("zlo zhi", line)
            if match:
                self.zlo = float(line.split()[0])
                self.zhi = float(line.split()[1])

                flag = 1
        f.close()

        subprocess.Popen(f"grep -A {self.num_atoms+1} Atoms {self.fname} | awk '{{print $1, $2, $3, $4, $5, $6}}' > atoms.txt", shell=True).wait()
        subprocess.Popen(f"grep -A {self.num_bonds+1} Bonds {self.fname} > bonds.txt", shell=True).wait()

        self.atoms = np.loadtxt("atoms.txt", skiprows=2)
        self.atoms = self.atoms[np.argsort(self.atoms[:,0])]

        self.bonds = np.loadtxt("bonds.txt", skiprows=2)
        self.bonds = self.bonds[np.argsort(self.bonds[:,2])]

        os.remove("atoms.txt")
        os.remove("bonds.txt")
        
    def probe(self):

        chain = self.atoms[self.atoms[:,1] == 1]
        types = np.unique(chain[:,2].astype(int))
        blocks = []
        for t in types:
            block = chain[chain[:,2] == t]
            blocks.append('-'.join([str(i) for i in block[:,2].astype(int)]))
        print('~'.join(blocks))

    def write_data(self, output='data_result.txt', header=None):

        f = open(os.path.join("results", output), "w")
        if header is not None:
            f.write(f"# {header}\n\n")
        else:
            f.write("# lmp data file\n\n")
        f.write('%d atoms\n' % self.num_atoms)
        f.write('%d atom types\n' % self.num_atomtypes)
        f.write('%d bonds\n' % self.num_bonds)
        f.write('%d bond types\n' % self.num_bondtypes)
        f.write('\n')
        f.write('%s %s xlo xhi\n' % (self.xlo, self.xhi))
        f.write('%s %s ylo yhi\n' % (self.ylo, self.yhi))
        f.write('%s %s zlo zhi\n' % (self.zlo, self.zhi))
        f.write('\n')
        f.write('Masses\n\n')
        for iind in range(self.num_atomtypes):
            f.write('%d 1\n' % (iind+1))
        f.write('\n')
        f.write('Atoms\n\n')
        for atom in self.atoms:
            f.write('%d\t%d\t%d\t%s\t%s\t%s\n' % (atom[0], atom[1], atom[2], atom[3], atom[4], atom[5]))
        f.write('\n')
        f.write('Bonds\n\n')
        for bond in self.bonds:
            f.write('%d\t%d\t%d\t%d\n' % (bond[0], bond[1], bond[2], bond[3]))
        f.close()

        if self.templates:
            f = open(os.path.join("results", self.script_name), "w")
            for line in self.templates:
                match = re.search("read_data", line)
                if match:
                    line = f"read_data\t\t{output}"
                f.write(line)
            f.close()

class Layer:
    def __init__(self, data, **kwargs):
        self.data = data
        self.args = kwargs
        
        if "film_types" in self.args:
            for ind, atom_type in enumerate(self.args["film_types"]):
                if ind == 0:
                    logic = self.data.atoms[:,2] == atom_type
                else:
                    logic = np.logical_or(logic, self.data.atoms[:,2] == atom_type)
            self.film = self.data.atoms[logic].copy()
            
        if "substrate_types" in self.args:
            for ind, substrate_type in enumerate(self.args["substrate_types"]):
                if ind == 0:
                    logic = self.data.atoms[:,2] == substrate_type
                else:
                    logic = np.logical_or(logic, self.data.atoms[:,2] == substrate_type)
            self.substrate = self.data.atoms[logic].copy()

        self.bonds = self.data.bonds.copy()

    def __truediv__(self, other):
        '''
        self on top of other; self over other; self / other
        '''
        data = LammpsData()
        data.xlo = self.data.xlo
        data.xhi = self.data.xhi
        data.ylo = self.data.ylo
        data.yhi = self.data.yhi
        data.zlo = self.data.zlo
        data.num_atoms = self.film.shape[0] + other.film.shape[0] + self.substrate.shape[0]
        data.num_atomtypes = len(self.args["film_types"]) + len(other.args["film_types"]) + len(self.args["substrate_types"])
        data.num_bonds = self.bonds.shape[0] + other.bonds.shape[0]
        data.num_bondtypes = len(np.unique(self.bonds[:,1])) + len(np.unique(other.bonds[:,1]))

        for old, new in zip(self.args["film_types"], self.args["res_film_types"]): 
            self.film[self.film[:,2] == old, 2] = new
        for old, new in zip(other.args["film_types"], other.args["res_film_types"]): 
            other.film[other.film[:,2] == old, 2] = new
            
        for old, new in zip(np.unique(self.bonds[:,1]), self.args["res_bond_types"]): 
            self.bonds[self.bonds[:,1] == old, 1] = new
        for old, new in zip(np.unique(other.bonds[:,1]), other.args["res_bond_types"]): 
            other.bonds[other.bonds[:,1] == old, 1] = new

        film_atom_types = self.args["res_film_types"] + other.args["res_film_types"]
        for old, new in zip(self.args["substrate_types"], np.array([1, 2]) + np.max(film_atom_types)):
            self.substrate[self.substrate[:,2] == old, 2] = new
        
        # self is put on top of other
        self.film[:,5] += np.max(other.film[:,5])

        # film: adjust atom-ID and atom-Type
        other.film[:,0] += np.max(self.film[:,0])
        other.film[:,1] += np.max(self.film[:,1])

        # substrate: adjust atom-ID and atom-Type
        self.substrate[:,0] = np.max(other.film[:,0]) + np.linspace(1, self.substrate.shape[0], self.substrate.shape[0])
        self.substrate[:,1] = np.max(other.film[:,1]) + 1

        # atoms: self.film, other film, and a substrate (self)
        data.atoms = np.row_stack((self.film, other.film, self.substrate))

        # bonds: adjust bond-ID and asscoiated atoms-ID's
        other.bonds[:,0] += np.max(self.bonds[:,0])
        other.bonds[:,2:] += self.film.shape[0]

        # bonds: self.bonds and other.bonds
        data.bonds = np.row_stack((self.bonds, other.bonds))

        # adjust zhi of a resulting film
        data.zhi = np.max(data.atoms[:,5]) + 20

        templates = []
        flag = 0
        f = open(os.path.join("templates", "t_welding.txt"), "r")
        while (flag == 0):
            line = f.readline()
            if len(line) > 0:
                templates.append(line)
            else:
                flag = 1
        f.close()

        data.templates = templates
        data.script_name = "in.welding.txt"

        return data

    def add_substrate(self, type_A=3, type_B=4, ratio=0.5):

        data = LammpsData()
        data.xlo = self.data.xlo
        data.xhi = self.data.xhi
        data.ylo = self.data.ylo
        data.yhi = self.data.yhi
        data.zlo = self.data.zlo
        data.zhi = self.data.zhi

        lx = self.data.xhi - self.data.xlo
        ly = self.data.yhi - self.data.ylo

        nx = int(lx/1)
        ny = int(math.ceil(ly/math.sqrt(3)*2))
        num_subsbeads = nx*ny
        coords = np.zeros([num_subsbeads, 6])
        rn = np.random.rand(len(coords))
        
        coords[:,0] = np.linspace(1, len(coords), len(coords)) + self.film.shape[0]
        for j in range(ny):
            for i in range(nx):
                if pow(-1,j) == 1:
                    coords[j*nx+i, 3:] = np.array([i, j*math.sqrt(3)/2, 0])
                else:
                    coords[j*nx+i, 3:] = np.array([0.5+i, j*math.sqrt(3)/2, 0])

        coords[:,1] = np.max(self.film[:,1]) + 1
        [coords[:,2], count_A, count_B] = self._det_type(type_A, type_B, rn)

        data.num_atoms = self.film.shape[0] + num_subsbeads
        data.num_atomtypes = len(self.args["film_types"]) + 2
        data.atoms = np.row_stack((self.film, coords))

        data.num_bonds = self.bonds.shape[0]
        data.num_bondtypes = len(np.unique(self.bonds[:,1]))
        data.bonds = self.bonds

        templates = []
        flag = 0
        f = open(os.path.join("templates", "t_glue.txt"), "r")
        while (flag == 0):
            line = f.readline()
            if len(line) > 0:
                templates.append(line)
            else:
                flag = 1
        f.close()

        data.templates = templates
        data.script_name = "in.glue.txt"

        return data

    def _det_type(self, type_A, type_B, rn):

        ttype = np.zeros(len(rn))
        count_A = 0
        count_B = 0
        for i in range(len(ttype)):
            ttype[i] = type_B
            count_B += 1
        else:
            ttype[i] = type_A
            count_A += 1

        return ttype, count_A, count_B

if __name__ == "__main__":
    main()
