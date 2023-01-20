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

    components = {
        0: {"N": 20, "f_A": 0.25, "n": 4800, "type_A": 1, "type_B": 2},
        #1: {"N": 22, "f_A": 0.5, "n": 4364, "type_A": 3, "type_B": 4},
        }

    # read in data (by self-avoiding random walk)
    #data = LammpsData("data_C20n2400_wo_substate.txt")

    # generate a layer
    #layer = Layer(data, film_types=[[1, 2],])

    # add a substrate
    #data_w_substrate = layer.add_substrate(substrate_types=[3, 4], ratio=0.5)

    # write the new data
    #data_w_substrate.write_data("data_tuned.txt", header="lmp data file")

    '''
    2. bilayer film
    '''
    # read in film(s) (already with their own substrates)
    data_C = LammpsData("data_C20n4800.txt")

    # generate layers
    C_top = Layer(data_C, film_types=[1, 2], substrate_types=[3, 4], res_film_types=[1, 2], res_bond_types=[1, 2, 3])
    C_bottom = Layer(data_C, film_types=[1, 2], substrate_types=[3, 4], res_film_types=[3, 4], res_bond_types=[4, 5, 6])

    # layer1 over layer2; use truediv to place one on top of the other
    ConC = C_top / C_bottom

    # write the new data
    ConC.write_data("data_tuned.txt", header="lmp data file of C20n4800 on C20n4800")


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
            for line in self.templates.format(output):
                f.write(line)
            f.close()

class Layer:
    def __init__(self, data, **kwargs):
        self.data = data
        self.args = kwargs
        
        if "film_types" in self.args:
            types = self.args["film_types"]
            if isinstance(types[0], list):
                types = [j for i in self.args["film_types"] for j in i]

            for ind, atom_type in enumerate(types):
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

        pairs = []
        types = [self.args["res_film_types"], other.args["res_film_types"]]
        print(types)
        for i, types_i in enumerate(types):
            type_Ai = types_i[0]
            type_Bi = types_i[1]
        
            pairs.append(f"pair_coeff\t\t\t{type_Ai} {type_Ai} lj/cut 1.0 1.0 2.5")
            pairs.append(f"pair_coeff\t\t\t{type_Ai} {type_Bi} lj/cut 1.0 1.0 2.5")
            pairs.append(f"pair_coeff\t\t\t{type_Bi} {type_Bi} lj/cut 1.0 1.0 2.5")

            for j, types_j in enumerate(types):
                type_Aj = types_j[0]
                type_Bj = types_j[1]

                if j > i:
                    pairs.append(f"pair_coeff\t\t\t{type_Ai} {type_Aj} lj/cut 1.0 1.0 2.5")
                    pairs.append(f"pair_coeff\t\t\t{type_Ai} {type_Bj} lj/cut 1.0 1.0 2.5")
                    pairs.append(f"pair_coeff\t\t\t{type_Bi} {type_Aj} lj/cut 1.0 1.0 2.5")
                    pairs.append(f"pair_coeff\t\t\t{type_Bi} {type_Bj} lj/cut 1.0 1.0 2.5")
            
        substrate_types = [np.max(types) + 1, np.max(types) + 2]
        for types in [self.args["res_film_types"], other.args["res_film_types"]]:
            type_A = types[0]
            type_B = types[1]
        
            for tp in substrate_types:
                pairs.append(f"pair_coeff\t\t\t{type_A} {tp} lj/cut 1.0 1.0 2.5")
                pairs.append(f"pair_coeff\t\t\t{type_B} {tp} lj/cut 1.0 1.0 2.5")

        pairs.append("\n")

        groups = []
        poly = []
        for types in [self.args["res_film_types"], other.args["res_film_types"]]:
            poly.append(str(types[0]))
            poly.append(str(types[1]))
        
        bottom = [str(tp) for tp in self.args["res_film_types"]]
        top = [str(tp) for tp in other.args["res_film_types"]]

        subs = [str(tp) for tp in substrate_types]
        
        groups.append(f"group\t\t\ttop type {' '.join(top)}")
        groups.append(f"group\t\t\tbottom type {' '.join(bottom)}")
        groups.append(f"group\t\t\tpoly type {' '.join(poly)}")
        groups.append(f"group\t\t\tsubs type {' '.join(subs)}")

        data.templates = t_w_1 + "\n".join(pairs) + "\n".join(groups) + t_w_2 
        data.script_name = "in.welding.txt"

        return data

    def add_substrate(self, substrate_types=[3, 4], ratio=0.5):

        
        #TODO: sDSA, grapho
        mode = None
        #

        data = LammpsData()
        data.xlo = self.data.xlo
        data.xhi = self.data.xhi
        data.ylo = self.data.ylo
        data.yhi = self.data.yhi
        data.zlo = self.data.zlo
        data.zhi = self.data.zhi

        self.substrate_types = substrate_types
        type_A = substrate_types[0]
        type_B = substrate_types[1]

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
        data.num_atomtypes = len([j for i in self.args["film_types"] for j in i]) + len(substrate_types)
        data.atoms = np.row_stack((self.film, coords))

        data.num_bonds = self.bonds.shape[0]
        data.num_bondtypes = len(np.unique(self.bonds[:,1]))
        data.bonds = self.bonds

        pairs = []
        for i, types_i in enumerate(self.args["film_types"]):
            type_Ai = types_i[0]
            type_Bi = types_i[1]
        
            pairs.append(f"pair_coeff\t\t\t{type_Ai} {type_Ai} lj/cut 1.0 1.0 2.5")
            pairs.append(f"pair_coeff\t\t\t{type_Ai} {type_Bi} lj/cut 1.0 1.0 2.5")
            pairs.append(f"pair_coeff\t\t\t{type_Bi} {type_Bi} lj/cut 1.0 1.0 2.5")

            for j, types_j in enumerate(self.args["film_types"]):
                type_Aj = types_j[0]
                type_Bj = types_j[1]

                if j > i:
                    pairs.append(f"pair_coeff\t\t\t{type_Ai} {type_Aj} lj/cut 1.0 1.0 2.5")
                    pairs.append(f"pair_coeff\t\t\t{type_Ai} {type_Bj} lj/cut 1.0 1.0 2.5")
                    pairs.append(f"pair_coeff\t\t\t{type_Bi} {type_Aj} lj/cut 1.0 1.0 2.5")
                    pairs.append(f"pair_coeff\t\t\t{type_Bi} {type_Bj} lj/cut 1.0 1.0 2.5")
            
        for types in self.args["film_types"]:
            type_A = types[0]
            type_B = types[1]
        
            for tp in self.substrate_types:
                pairs.append(f"pair_coeff\t\t\t{type_A} {tp} lj/cut 1.0 1.0 2.5")
                pairs.append(f"pair_coeff\t\t\t{type_B} {tp} lj/cut 1.0 1.0 2.5")

        pairs.append("\n")

        pairs_mod = []
        pairs_mod.append("\n")
        for types in self.args["film_types"]:
            type_A = types[0]
            type_B = types[1]
        
            for tp in self.substrate_types:
                pairs_mod.append(f"pair_coeff\t\t\t{type_A} {tp} lj/cut 0.745 1.0 2.5")
                pairs_mod.append(f"pair_coeff\t\t\t{type_B} {tp} lj/cut 0.745 1.0 2.5")

        pairs_mod.append("\n")

        groups = []
        poly = []
        for types in self.args["film_types"]:
            poly.append(str(types[0]))
            poly.append(str(types[1]))
        
        subs = [str(tp) for tp in self.substrate_types]
        
        groups.append(f"group\t\t\tpoly type {' '.join(poly)}")
        groups.append(f"group\t\t\tsubs type {' '.join(subs)}")

        data.templates = t_g_1 + "\n".join(pairs) + "\n".join(groups) + t_g_2 + "\n".join(pairs_mod) + t_g_3
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

t_g_1 = '''
# lammps input file
# Variable
variable		T_start equal 1.2
variable		T_end equal 1.2
variable		T_damp equal 10.0

# Initialization
units			lj
boundary		p p f
atom_style		bond

# Forcefield and coordinates data
neighbor		0.3 bin
neigh_modify	delay 0 one 2000 page 20000
bond_style		fene
pair_style		hybrid lj/cut 2.5
pair_modify		shift yes 

read_data		{}

bond_coeff		* 30.0 1.5 1.0 1.0

special_bonds	fene angle no dihedral no lj/coul 0 1 1

comm_style		tiled

pair_coeff      * * none
'''

t_g_2 = '''
#####################################################
fix             zwall_lo poly wall/reflect zlo EDGE
fix             zwall_hi poly wall/reflect zhi EDGE
fix				zeroforce subs setforce 0 0 0

minimize		1.0e-10 1.0e-10 100000 10000000

write_data		data_min
write_restart	restart_min

velocity		all create ${{T_start}} 900531 dist gaussian

reset_timestep  0
'''

t_g_3 = '''
timestep		0.006

fix				bal all balance 100000 1.1 rcb out balancing
fix				zeromomentum poly momentum 100 linear 1 1 1

dump			simple all custom 10000 dump.* id mol type x y z 

fix             1 poly nvt temp ${{T_start}} ${{T_end}} ${{T_damp}}

thermo_style	custom step temp press ke pe epair ebond pxx pyy pzz vol density
thermo          100 

run				1000000

unfix			1
unfix			zeromomentum

write_data		data_preequil
write_restart	restart_preequil

print "All done"
'''

t_w_1 = '''
# lammps input file 
# Variable
variable		T_start equal 1.2
variable		T_end equal 1.2
variable		T_damp equal 10.0

# Initialization
units			lj
boundary		p p f
atom_style		bond

# Forcefield and coordinates data
neighbor		0.3 bin
neigh_modify	delay 0 one 2000 page 20000
bond_style		fene
pair_style		hybrid lj/cut 2.5
pair_modify		shift yes 

read_data		{}

bond_coeff		* 30.0 1.5 1.0 1.0

special_bonds	fene angle no dihedral no lj/coul 0 1 1

comm_style		tiled

pair_coeff      * * none
'''

t_w_2 = '''
#####################################################
reset_timestep  0

velocity		all create ${{T_start}} 900531 dist gaussian

timestep		0.006

fix				zwall_lo poly wall/reflect zlo EDGE
fix				zwall_hi poly wall/reflect zhi EDGE
fix				bal all balance 100000 1.1 rcb out balancing
fix				zeroforce subs setforce 0 0 0
fix				zeroforce bottom setforce 0 0 0
fix				zeromomentum poly momentum 100 linear 1 1 1

dump			simple all custom 10000 dump.* id mol type x y z 

fix				1 top nvt temp ${{T_start}} ${{T_end}} ${{T_damp}}

thermo_style	custom step temp press ke pe epair ebond pxx pyy pzz vol density
thermo			100 

run				100000

unfix			1
unfix			zeromomentum

write_data		data_preequil
write_restart	restart_preequil

print "All done"
'''

if __name__ == "__main__":
    main()
