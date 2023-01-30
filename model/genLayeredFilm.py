#!/usr/bin/env python3

import numpy as np
import re
import subprocess
import os
import math


class LammpsData:

    '''Lammps data'''
    def __init__(self, fname=None):
        if fname is not None:
            self.fname = os.path.join("ingredients", fname)
            self.read_data()
            self.has_substrate = False

            self.probe()
            self.templates = None

    def read_data(self):
        '''
        read in atom coords and bond topology information
        '''

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

        components = {}
        substrate = {}

        trj = self.atoms

        mol_ids = np.unique(trj[:,1]).astype(int)
        self.mol_ids = mol_ids
        mols = np.empty(len(mol_ids))
        self.mols = mols

        for mol_ind, mol_id in enumerate(mol_ids):
            mol = trj[trj[:,1] == mol_id]
            types = np.unique(mol[:,2]).astype(int)

            N = mol.shape[0]
            N_A = mol[mol[:,2] == types[0]].shape[0]
            f_A = N_A/N

            if len(components) == 0:
                ind = len(components)
                components.update({
                    ind: {'N': N, 'f_A': f_A, "block_types": {'A': types[0], 'B': types[1]}, 'n': 0}
                    })
                mols[mol_ind] = ind

            flag = 0
            for key in components:
                component = components.get(key)
                if (component['N'] == N and component['f_A'] == f_A and component["block_types"].get('A') == types[0] and component["block_types"].get('B') == types[1]):
                    mols[mol_ind] = key
                    component.update({'n': component.get('n') + 1})
                    flag = 1
                    break

            if flag == 0:
                ind = len(components)
                components.update({
                    ind: {'N': N, 'f_A': f_A, "block_types": {'A': types[0], 'B': types[1]}, 'n': 1}
                    })
                mols[mol_ind] = ind

        for key in components:
            component = components.get(key)
            blocks = []
            
            if component['N'] < 100:

                blocks.append('-'.join(['A' for _ in range(int(component['N']*component['f_A']))]))
                blocks.append('-'.join(['B' for _ in range(int(component['N']*(1-component['f_A'])))]))
            
                print(f"** component {key}")
                print(f"* {component}")
                #print(f"* structure: {'~'.join(blocks)}")

            else:
                key_substrate = key
                self.has_substrate = True

        
        for key in components:
            component = components.get(key)

            type_A = component.get('block_types').get('A')
            type_B = component.get('block_types').get('B')

            # get a representative mol
            mol_id = min(mol_ids[mols == key])

            mol = trj[trj[:,1] == mol_id]

            if len(mol) < 100:

                bond_types = {}

                for atom_ind in range(mol.shape[0]-1):
                    atom_1 = mol[atom_ind, 0]
                    atom_2 = mol[atom_ind + 1, 0]
                
                    logic = np.logical_and(self.bonds[:,2] == atom_1, self.bonds[:,3] == atom_2)

                    if mol[atom_ind, 2] == type_A and mol[atom_ind + 1, 2] == type_A:
                        bond_types.update({'AA': self.bonds[logic, 1][0].astype(int)})

                    if mol[atom_ind, 2] == type_A and mol[atom_ind + 1, 2] == type_B:
                        bond_types.update({'AB': self.bonds[logic, 1][0].astype(int)})

                    if mol[atom_ind, 2] == type_B and mol[atom_ind + 1, 2] == type_B:
                        bond_types.update({'BB': self.bonds[logic, 1][0].astype(int)})

                component.update({'bond_types': bond_types})

        if self.has_substrate:
            substrate.update(components.get(key_substrate))
            components.pop(key_substrate)

        self.components = components
        self.substrate = substrate


    def tailor(self, component_key, **kwargs):
        
        component = self.components.get(component_key)

        if 'block_types' in kwargs:
            type_A = kwargs['block_types'].get('A')
            type_B = kwargs['block_types'].get('B')
        else:
            type_A = component.get('block_types').get('A')
            type_B = component.get('block_types').get('B')

        if 'bond_types' in kwargs:
            type_AA = kwargs['bond_types'].get('AA')
            type_AB = kwargs['bond_types'].get('AB')
            type_BB = kwargs['bond_types'].get('BB')
        else:
            type_AA = component.get('bond_types').get('AA')
            type_AB = component.get('bond_types').get('AA')
            type_BB = component.get('bond_types').get('BB')

        if 'f_A' in kwargs:
            f_A = kwargs['f_A']
        else:
            f_A = component.get('f_A')

        N = component.get('N')

        self._atoms = self.atoms.copy()
        self._bonds = self.bonds.copy()
        if ('block_types' in kwargs) or ('bond_types' in kwargs) or ('f_A' in kwargs):

            mol_ids = self.mol_ids[self.mols == component_key]

            temp_types = np.concatenate((np.ones(int(N*f_A))*type_A, np.ones(int(N*(1-f_A)))*type_B))

            for mol_ind, mol_id in enumerate(mol_ids):
                self._atoms[self.atoms[:,1] == mol_id, 2] = temp_types

                mol = self._atoms[self.atoms[:,1] == mol_id]

                for atom_ind in range(mol.shape[0]-1):
                    atom_1 = mol[atom_ind, 0]
                    atom_2 = mol[atom_ind + 1, 0]
                
                    # method 1)
                    #cond_1 = np.argwhere(self.bonds[:,2] == atom_1)
                    #cond_2 = np.argwhere(self.bonds[:,3] == atom_2)
                    #bond_ind = np.intersect1d(cond_1, cond_2).astype(int)[0]

                    # method 2)
                    bond_ind = np.argwhere(self.bonds[:,2] == atom_1).astype(int)[0]
                    
                    # method 3)
                    #bond_ind = (N-1) * mol_ind + atom_ind
                    #print(bond_ind)

                    if mol[atom_ind, 2] == type_A and mol[atom_ind + 1, 2] == type_A:
                        self._bonds[bond_ind, 1] = type_AA

                    if mol[atom_ind, 2] == type_A and mol[atom_ind + 1, 2] == type_B:
                        self._bonds[bond_ind, 1] = type_AB

                    if mol[atom_ind, 2] == type_B and mol[atom_ind + 1, 2] == type_B:
                        self._bonds[bond_ind, 1] = type_BB
                    
        self.atoms = self._atoms
        self.bonds = self._bonds
        self.components[component_key].update({'f_A': f_A})
        self.components[component_key].update({'block_types': {'A': type_A, 'B': type_B}})
        self.components[component_key].update({'bond_types': {'AA': type_AA, 'AB': type_AB, 'BB': type_BB}})
        self.num_atomtypes = len(np.unique(self.atoms[:,2]))
        self.num_bondtypes = len(np.unique(self.bonds[:,1]))



    def write(self, output='data_result.txt', header=None):

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
        
        if self.templates is not None:
            f = open(os.path.join("results", self.script_name), "w")
            for line in self.templates.format(output):
                f.write(line)
            f.close()

class Layer:
    def __init__(self, data, blend=None, **kwargs):
        self.data = data
        self.args = kwargs
            

        # instead of starting with an empty dict, if given [1, 2], insert 1 and 2 into a list, pop 2 and update 1
        if blend is not None:
            components = {}

            sub_components = []

            for component_key in self.data.components:
                
                if component_key in blend:
                    sub_components.append(self.data.components.get(component_key))
                else:
                    components.update({component_key: self.data.components.get(component_key)})

            components.update({blend[0]: sub_components})
            self.data.components = components


        types = []
        bond_types = []
        substrate_types = {}

        for c_ind, component_key in enumerate(self.data.components):
            component = self.data.components.get(component_key)

            if isinstance(component, list):
                # blend
                for sub_component in component:
                    types.append(sub_component.get('block_types'))
                    bond_types.append(sub_component.get('bond_types'))

            elif isinstance(component, dict):
                # pure
                # a layer could be a layered film
                types.append(component.get('block_types'))
                bond_types.append(component.get('bond_types'))

        if self.data.substrate:
            # substrate

            substrate_types.update(self.data.substrate.get('block_types'))

        self.types = types
        self.bond_types = bond_types

        for ind, atom_type in enumerate(types):
            if ind == 0:
                logic_A = self.data.atoms[:,2] == atom_type.get('A')
                logic_B = self.data.atoms[:,2] == atom_type.get('B')

                logic = np.logical_or(logic_A, logic_B)

            else:
                logic_A = self.data.atoms[:,2] == atom_type.get('A')
                logic_B = self.data.atoms[:,2] == atom_type.get('B')

                logic = np.logical_or(logic, np.logical_or(logic_A, logic_B))

        self.film = self.data.atoms[logic].copy()
            
        if self.data.substrate:

            self.substrate_types = substrate_types

            logic_A = self.data.atoms[:,2] == self.data.substrate.get('block_types').get('A')
            logic_B = self.data.atoms[:,2] == self.data.substrate.get('block_types').get('B')

            logic = np.logical_or(logic_A, logic_B)

            self.substrate = self.data.atoms[logic].copy()

        self.bonds = self.data.bonds.copy()


    def __truediv__(self, other):
        '''
        self on top of other; self over other; self / other
        '''
        # TODO: change the return type to Layer

        data = LammpsData()
        data.xlo = self.data.xlo
        data.xhi = self.data.xhi
        data.ylo = self.data.ylo
        data.yhi = self.data.yhi
        data.zlo = self.data.zlo
        data.num_atoms = self.film.shape[0] + other.film.shape[0] + self.substrate.shape[0]
        data.num_bonds = self.bonds.shape[0] + other.bonds.shape[0]

        data.num_atomtypes = len(self.types)*2 + len(other.types)*2 + len(other.substrate_types)
        data.num_bondtypes = len(self.bond_types)*3 + len(other.bond_types)*3

        data.components = other.data.components
        offset = 0
        offset_bond = 0
        self._film = self.film.copy()
        self._bonds = self.bonds.copy()
        for component_key in self.data.components:
            component = self.data.components.get(component_key)

            if isinstance(component, list):
                # blend
                for sub_component in component:
                    old_A = sub_component.get('block_types').get('A')
                    old_B = sub_component.get('block_types').get('B')
                    new_A = len(other.types)*2 + offset + 1
                    new_B = len(other.types)*2 + offset + 2
                    sub_component.update({'block_types': {'A': new_A, 'B': new_B}})
                    
                    self._film[self.film[:,2] == old_A, 2] = new_A
                    self._film[self.film[:,2] == old_B, 2] = new_B

                    old_AA = sub_component.get('bond_types').get('AA')
                    old_AB = sub_component.get('bond_types').get('AB')
                    old_BB = sub_component.get('bond_types').get('BB')
                    new_AA = len(other.bond_types)*3 + offset_bond + 1 
                    new_AB = len(other.bond_types)*3 + offset_bond + 3
                    new_BB = len(other.bond_types)*3 + offset_bond + 2

                    sub_component.update({'bond_types': {'AA': new_AA, 'AB': new_AB, 'BB': new_BB}})

                    self._bonds[self.bonds[:,1] == old_AA, 1] = new_AA
                    self._bonds[self.bonds[:,1] == old_AB, 1] = new_AB
                    self._bonds[self.bonds[:,1] == old_BB, 1] = new_BB

                    offset += 2
                    offset_bond += 3

            elif isinstance(component, dict):
                # pure
                # a layer could be a layered film
                print('pure')
                old_A = component.get('block_types').get('A')
                old_B = component.get('block_types').get('B')
                new_A = len(other.types)*2 + offset + 1
                new_B = len(other.types)*2 + offset + 2

                component.update({'block_types': {'A': new_A, 'B': new_B}})
                    
                self._film[self.film[:,2] == old_A, 2] = new_A
                self._film[self.film[:,2] == old_B, 2] = new_B

                old_AA = component.get('bond_types').get('AA')
                old_AB = component.get('bond_types').get('AB')
                old_BB = component.get('bond_types').get('BB')
                new_AA = len(other.bond_types)*3 + offset_bond + 1 
                new_AB = len(other.bond_types)*3 + offset_bond + 3
                new_BB = len(other.bond_types)*3 + offset_bond + 2

                component.update({'bond_types': {'AA': new_AA, 'AB': new_AB, 'BB': new_BB}})

                self._bonds[self.bonds[:,1] == old_AA, 1] = new_AA
                self._bonds[self.bonds[:,1] == old_AB, 1] = new_AB
                self._bonds[self.bonds[:,1] == old_BB, 1] = new_BB

                offset += 2
                offset_bond += 3

            data.components.update({len(data.components): component})
            self.data.components.update({component_key: component})
                
        old_A = other.substrate_types.get('A')
        new_A = len(other.types)*2 + offset + 1
        old_B = other.substrate_types.get('B')
        new_B = len(other.types)*2 + offset + 2
        other.data.substrate.update({'block_types': {'A': new_A, 'B': new_B}})
        other.substrate[other.substrate[:,2] == old_A, 2] = new_A
        other.substrate[other.substrate[:,2] == old_B, 2] = new_B
        data.substrate = other.data.substrate

        self._film[:,0] += max(other.film[:,0])
        self._film[:,1] += max(other.film[:,1])
        self._film[:,5] += max(other.film[:,5])
        other.substrate[:,0] = np.linspace(1, other.substrate.shape[0], other.substrate.shape[0]).astype(int) + max(self._film[:,0])
        other.substrate[:,1] = max(self._film[:,1]) + 1

        data.atoms = np.row_stack((other.film, self._film, other.substrate))

        self._bonds[:,0] += max(other.bonds[:,0])
        self._bonds[:,2:] += max(other.film[:,0])
        data.bonds = np.row_stack((other.bonds, self._bonds))

        layer = Layer(data)

        # adjust zhi of a resulting film
        data.zhi = np.max(data.atoms[:,5]) + 20

        pairs = []
        for i, types_i in enumerate(layer.types):
            type_Ai = types_i.get('A')
            type_Bi = types_i.get('B')
        
            pairs.append(f"pair_coeff\t\t\t{type_Ai} {type_Ai} lj/cut 1.0 1.0 2.5")
            pairs.append(f"pair_coeff\t\t\t{type_Ai} {type_Bi} lj/cut 1.0 1.0 2.5")
            pairs.append(f"pair_coeff\t\t\t{type_Bi} {type_Bi} lj/cut 1.0 1.0 2.5")

            for j, types_j in enumerate(layer.types):
                type_Aj = types_j.get('A')
                type_Bj = types_j.get('B')

                if j > i:
                    pairs.append(f"pair_coeff\t\t\t{type_Ai} {type_Aj} lj/cut 1.0 1.0 2.5")
                    pairs.append(f"pair_coeff\t\t\t{type_Ai} {type_Bj} lj/cut 1.0 1.0 2.5")
                    pairs.append(f"pair_coeff\t\t\t{type_Bi} {type_Aj} lj/cut 1.0 1.0 2.5")
                    pairs.append(f"pair_coeff\t\t\t{type_Bi} {type_Bj} lj/cut 1.0 1.0 2.5")
            
        for types in layer.types:
            type_A = types.get('A')
            type_B = types.get('B')
        
            tp_A = layer.substrate_types.get('A')
            tp_B = layer.substrate_types.get('B')
            pairs.append(f"pair_coeff\t\t\t{type_A} {tp_A} lj/cut 0.745 1.0 2.5")
            pairs.append(f"pair_coeff\t\t\t{type_B} {tp_A} lj/cut 0.745 1.0 2.5")
            pairs.append(f"pair_coeff\t\t\t{type_A} {tp_B} lj/cut 0.745 1.0 2.5")
            pairs.append(f"pair_coeff\t\t\t{type_B} {tp_B} lj/cut 0.745 1.0 2.5")
        
        pairs.append("\n")

        groups = []
        poly = []
        top = []
        bottom = []
        for types in layer.types:
            poly.append(str(types.get('A')))
            poly.append(str(types.get('B')))
        
        for types in other.types:
            bottom.append(str(types.get('A')))
            bottom.append(str(types.get('B')))

        for component_key in self.data.components:
            component = self.data.components.get(component_key)
            if isinstance(component, list):
                for sub_component in component:
                    top.append(str(sub_component.get('block_types').get('A')))
                    top.append(str(sub_component.get('block_types').get('B')))

            else:
                top.append(str(component.get('block_types').get('A')))
                top.append(str(component.get('block_types').get('B')))

        subs = [str(tp) for tp in layer.substrate_types.values()]
        
        groups.append(f"group\t\t\tpoly type {' '.join(poly)}")
        groups.append(f"group\t\t\tsubs type {' '.join(subs)}")
 
        groups.append(f"group\t\t\ttop type {' '.join(top)}")
        groups.append(f"group\t\t\tbottom type {' '.join(bottom)}")
        groups.append(f"group\t\t\tpoly type {' '.join(poly)}")
        groups.append(f"group\t\t\tsubs type {' '.join(subs)}")

        layer.data.templates = t_w_1 + "\n".join(pairs) + "\n".join(groups) + t_w_2 
        layer.data.script_name = "in.welding.txt"
        layer.data.has_substrate = True

        return layer

    def add_substrate(self, substrate_types={'A': 3, 'B': 4}, ratio=0.5):

        #TODO: sDSA, grapho
        mode = None
        #

        self.substrate_types = substrate_types
        self.substrate_ratio = ratio
        type_A = substrate_types.get('A')
        type_B = substrate_types.get('B')

        lx = self.data.xhi - self.data.xlo
        ly = self.data.yhi - self.data.ylo

        nx = int(lx/1)
        ny = int(math.ceil(ly/math.sqrt(3)*2))
        num_subsbeads = nx*ny
        substrate = np.zeros([num_subsbeads, 6])
        self.num_subsbeads = num_subsbeads
                
        substrate[:,0] = np.linspace(1, len(substrate), len(substrate)) + self.film.shape[0]
        for j in range(ny):
            for i in range(nx):
                if pow(-1,j) == 1:
                    substrate[j*nx+i, 3:] = np.array([i, j*math.sqrt(3)/2, 0])
                else:
                    substrate[j*nx+i, 3:] = np.array([0.5+i, j*math.sqrt(3)/2, 0])

        substrate[:,1] = np.max(self.film[:,1]) + 1
        [substrate[:,2], count_A, count_B] = self._det_type()

        self.data.num_atoms = self.film.shape[0] + num_subsbeads
        self.data.num_atomtypes = len(self.types)*2 + len(substrate_types)
        self.data.atoms = np.row_stack((self.film, substrate))
        self.substrate = substrate

        self.data.num_bonds = self.bonds.shape[0]
        self.data.num_bondtypes = len(np.unique(self.bonds[:,1]))
        self.data.bonds = self.bonds

        pairs = []
        for i, types_i in enumerate(self.types):
            type_Ai = types_i.get('A')
            type_Bi = types_i.get('B')
        
            pairs.append(f"pair_coeff\t\t\t{type_Ai} {type_Ai} lj/cut 1.0 1.0 2.5")
            pairs.append(f"pair_coeff\t\t\t{type_Ai} {type_Bi} lj/cut 1.0 1.0 2.5")
            pairs.append(f"pair_coeff\t\t\t{type_Bi} {type_Bi} lj/cut 1.0 1.0 2.5")

            for j, types_j in enumerate(self.types):
                type_Aj = types_j.get('A')
                type_Bj = types_j.get('B')

                if j > i:
                    pairs.append(f"pair_coeff\t\t\t{type_Ai} {type_Aj} lj/cut 1.0 1.0 2.5")
                    pairs.append(f"pair_coeff\t\t\t{type_Ai} {type_Bj} lj/cut 1.0 1.0 2.5")
                    pairs.append(f"pair_coeff\t\t\t{type_Bi} {type_Aj} lj/cut 1.0 1.0 2.5")
                    pairs.append(f"pair_coeff\t\t\t{type_Bi} {type_Bj} lj/cut 1.0 1.0 2.5")
            
        for types in self.types:
            type_A = types.get('A')
            type_B = types.get('B')
        
            tp_A = substrate_types.get('A')
            tp_B = substrate_types.get('B')
            pairs.append(f"pair_coeff\t\t\t{type_A} {tp_A} lj/cut 1.0 1.0 2.5")
            pairs.append(f"pair_coeff\t\t\t{type_B} {tp_A} lj/cut 1.0 1.0 2.5")
            pairs.append(f"pair_coeff\t\t\t{type_A} {tp_B} lj/cut 1.0 1.0 2.5")
            pairs.append(f"pair_coeff\t\t\t{type_B} {tp_B} lj/cut 1.0 1.0 2.5")

        pairs.append("\n")

        pairs_mod = []
        pairs_mod.append("\n")
        for types in self.types:
            type_A = types.get('A')
            type_B = types.get('B')
        
            tp_A = substrate_types.get('A')
            tp_B = substrate_types.get('B')
            pairs_mod.append(f"pair_coeff\t\t\t{type_A} {tp_A} lj/cut 0.745 1.0 2.5")
            pairs_mod.append(f"pair_coeff\t\t\t{type_B} {tp_A} lj/cut 0.745 1.0 2.5")
            pairs_mod.append(f"pair_coeff\t\t\t{type_A} {tp_B} lj/cut 0.745 1.0 2.5")
            pairs_mod.append(f"pair_coeff\t\t\t{type_B} {tp_B} lj/cut 0.745 1.0 2.5")

        pairs_mod.append("\n")

        groups = []
        poly = []
        for types in self.types:
            poly.append(str(types.get('A')))
            poly.append(str(types.get('B')))
        
        subs = [str(tp) for tp in substrate_types.values()]
        
        groups.append(f"group\t\t\tpoly type {' '.join(poly)}")
        groups.append(f"group\t\t\tsubs type {' '.join(subs)}")

        self.data.templates = t_g_1 + "\n".join(pairs) + "\n".join(groups) + t_g_2 + "\n".join(pairs_mod) + t_g_3
        self.data.script_name = "in.glue.txt"

        self.data.has_substrate = True

    def _det_type(self):

        rn = np.random.rand(self.num_subsbeads)
        ttype = np.zeros(len(rn))
        
        B = rn > self.substrate_ratio
        A = rn <= self.substrate_ratio
        ttype[B] = self.substrate_types.get('B')
        ttype[A] = self.substrate_types.get('A')
        
        return ttype, sum(A), sum(B)

t_g_1 = '''# lammps input file
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
