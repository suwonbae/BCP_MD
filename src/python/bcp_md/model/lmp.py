#!/usr/bin/env python3

import numpy as np
import re
import subprocess
import os
import math

from . import rw
from .rw import *

__all__ = ['Data', 'Layer', 'Script']

class Data:
    """
    Attributes:
        fname (str): lmp data file name
    """
    def __init__(self, fname=None):
        """
        Data constructor
        """

        self.has_substrate = False
        self.templates = None

        if fname is not None:
            #self.fname = os.path.join("ingredients", fname)
            self.fname = fname
            self.read_data()
            self.probe()


    def read_data(self):
        """
        read in atom coords and bond topology information
        """

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
        self._atoms = self.atoms.copy()

        self.bonds = np.loadtxt("bonds.txt", skiprows=2)
        self.bonds = self.bonds[np.argsort(self.bonds[:,2])]
        self._bonds = self.bonds.copy()

        os.remove("atoms.txt")
        os.remove("bonds.txt")


    def generate(self, components, sim_box, bc):
        """
        generate a self-avoiding random walk structure
        Parameters:
        ----------
        components (dict): components to be used to construct layers and films
        sim_box (list): list of lower and upper bounds in each axis; [[xlo, xhi], [ylo, yhi], [zlo, zhi]]
        bc (list): list of boundary conditions in each axis; ['p', 'p', 'f']

        Returns:
        -------
        N/A
        """

        params = []

        params.append(f"components {len(components)}")

        for key in components:
            component = components.get(key)

            params.append(f"N f_A A B n {component.get('N')} {component.get('f_A')} {component.get('block_types').get('A')} {component.get('block_types').get('B')} {component.get('n')}")

        if len(sim_box) != 3: raise Exception("sim_box has to be 3-dim")
        params.append(f"xlo xhi {sim_box[0][0]} {sim_box[0][1]}")
        params.append(f"ylo yhi {sim_box[1][0]} {sim_box[1][1]}")
        params.append(f"zlo zhi {sim_box[2][0]} {sim_box[2][1]}")

        if len(bc) != 3: raise Exception("boundary condition has to be 3-dim")
        params.append(f"boundary {bc[0]} {bc[1]} {bc[2]}")

        if (components is not None) and (sim_box is not None) and (bc is not None):
            print("random walk")

            f = open('in.parameters', 'w')
            f.write('\n'.join(params))
            f.close()

            # call a function using pybind11
            fn_rw()

            self.fname = 'data.txt'
            self.read_data()
            self.probe()


    def probe(self):
        """
        probe the given lmp data file and get components
        """

        components = {}

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
            components.update({'substrate': components.get(key_substrate)})
            components.pop(key_substrate)

        self.components = components


    def tailor(self, component_key, **kwargs):
        """
        tailor f_A, block_types, or/and bond_types for a component
        Parameters:
        ----------
        component_key (int/string): component key that refers to the component to be tailored
        kwargs: only 'f_A', 'block_types', and 'bond_types' work
            block_types (dict): new types of constituent beads; {'A': 1, 'B': 2}
            bond_types (dict): new types of bonds; {'AA': 1, 'AB': 3, 'BB': 2}
            f_A (float): new fraction of minority A; architecture of linear BCP

        Returns:
        -------
        N/A
        """

        component = self.components.get(component_key)

        if component_key == 'substrate':
            if 'block_types' in kwargs:
                old_A = component.get('block_types').get('A')
                new_A = kwargs['block_types'].get('A')
                old_B = component.get('block_types').get('B')
                new_B = kwargs['block_types'].get('B')

                self.atoms[self._atoms[:,2] == old_A, 2] = new_A
                self.atoms[self._atoms[:,2] == old_B, 2] = new_B

            self.components[component_key].update({'block_types': {'A': new_A, 'B': new_B}})
            self.num_atomtypes = len(np.unique(self.atoms[:,2]))
            self.num_bondtypes = len(np.unique(self.bonds[:,1]))

        else:
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
    
            
            if ('block_types' in kwargs) or ('bond_types' in kwargs) or ('f_A' in kwargs):

                mol_ids = self.mol_ids[self.mols == component_key]

                temp_types = np.concatenate((np.ones(int(N*f_A))*type_A, np.ones(int(N*(1-f_A)))*type_B))

                for mol_ind, mol_id in enumerate(mol_ids):
                    self.atoms[self._atoms[:,1] == mol_id, 2] = temp_types

                    mol = self.atoms[self._atoms[:,1] == mol_id]

                    for atom_ind in range(mol.shape[0]-1):
                        atom_1 = mol[atom_ind, 0]
                        atom_2 = mol[atom_ind + 1, 0]
                
                        # method 1)
                        #cond_1 = np.argwhere(self._bonds[:,2] == atom_1)
                        #cond_2 = np.argwhere(self._bonds[:,3] == atom_2)
                        #bond_ind = np.intersect1d(cond_1, cond_2).astype(int)[0]

                        # method 2)
                        bond_ind = np.argwhere(self._bonds[:,2] == atom_1).astype(int)[0]
                    
                        # method 3)
                        #bond_ind = (N-1) * mol_ind + atom_ind

                        if mol[atom_ind, 2] == type_A and mol[atom_ind + 1, 2] == type_A:
                            self.bonds[bond_ind, 1] = type_AA

                        if mol[atom_ind, 2] == type_A and mol[atom_ind + 1, 2] == type_B:
                            self.bonds[bond_ind, 1] = type_AB

                        if mol[atom_ind, 2] == type_B and mol[atom_ind + 1, 2] == type_B:
                            self.bonds[bond_ind, 1] = type_BB
                    
            self.components[component_key].update({'f_A': f_A})
            self.components[component_key].update({'block_types': {'A': type_A, 'B': type_B}})
            self.components[component_key].update({'bond_types': {'AA': type_AA, 'AB': type_AB, 'BB': type_BB}})
            self.num_atomtypes = len(np.unique(self.atoms[:,2]))
            self.num_bondtypes = len(np.unique(self.bonds[:,1]))


    def modify_component(self, from_c, to_c):
        """
        modify component_key since the layer ordering may not correctly assigned when reading in an lmp data file
        Parameters:
        ----------
        from_c (list): list of current ordering; [0, 1]
        to_c (list): list of new ordering; [1, 0]

        Returns:
        -------
        N/A
        """

        components = self.components.copy()

        components_tmp = {}
            
        old_keys = from_c # ex: [0, 1]
        new_keys = to_c #ex: [1, 0]

        for old_key, new_key in zip(old_keys, new_keys):
            print(old_key, new_key)
            components_tmp.update({new_key: components.get(old_key)})
            
        components_tmp.update({'substrate': components.get('substrate')})
        print("tmp", components_tmp)
        components = components_tmp

        self.components = components


    def write(self, output='data_result.txt', header=None):
        """
        write an lmp data file
        Parameters:
        ----------
        output (str): name of lmp data file to be written
        header (str): header/comment to be written on the top  

        Returns:
        -------
        N/A
        """

        #f = open(os.path.join("results", output), "w")
        f = open(output, "w")
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
            #f = open(os.path.join("results", self.script_name), "w")
            f = open(self.script_name, "w")

            for line in self.templates.format(output):
                f.write(line)
            f.close()


class Layer:
    """
    Attributes:
        data (Data): an instance of the lmp data type
        blend (list): list of components that are blended
    """
    def __init__(self, data, blend=None, **kwargs):
        """
        Layer constructor
        """

        if not isinstance(data, Data): raise Exception("data has to be an lmp data")

        self.data = data
        self.args = kwargs
        
        # TODO:
        # as is okay, but instead of starting with an empty dict
        # if given [1, 2], insert 1 and 2 into a list, pop 2 and update 1
        if blend is not None:
            components = {}

            sub_components = []

            for component_key, component in self.data.components.items():
                
                if component_key in blend:
                    sub_components.append(component)
                else:
                    components.update({component_key: component})

            components.update({blend[0]: sub_components})
            self.data.components = components

        types = []
        bond_types = []
        substrate_types = {}

        for component_key in self.data.components:
            component = self.data.components.get(component_key)

            if component_key != 'substrate':
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

        if self.data.has_substrate:
            # substrate

            substrate_types.update(self.data.components.get('substrate').get('block_types'))

        self.types = types
        self.bond_types = bond_types
        print('hi', types)

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
            
        if self.data.has_substrate:

            self.substrate_types = substrate_types

            logic_A = self.data.atoms[:,2] == self.data.components.get('substrate').get('block_types').get('A')
            logic_B = self.data.atoms[:,2] == self.data.components.get('substrate').get('block_types').get('B')

            logic = np.logical_or(logic_A, logic_B)

            self.substrate = self.data.atoms[logic].copy()

        self.bonds = self.data.bonds.copy()


    def __truediv__(self, other):
        """
        generate layered structure by putting self on top of other; self over other; self / other

        this is a defined operator, not a function
        both operands should be of the Layer type
        Parameters:
        ----------
        N/A
        
        Returns:
        -------
        layer (Layer): a new instance of the Layer type
        """

        if not isinstance(self, Layer): raise Exception("the numerator has to be a layer")
        if not isinstance(other, Layer): raise Exception("the denominator has to be a layer")

        data = Data()
        data.xlo = self.data.xlo
        data.xhi = self.data.xhi
        data.ylo = self.data.ylo
        data.yhi = self.data.yhi
        data.zlo = self.data.zlo
        data.num_atoms = self.film.shape[0] + other.film.shape[0] + self.substrate.shape[0]
        data.num_bonds = self.bonds.shape[0] + other.bonds.shape[0]

        data.num_atomtypes = len(self.types)*2 + len(other.types)*2 + len(other.substrate_types)
        print(self.types, other.types, other.substrate_types)
        data.num_bondtypes = len(self.bond_types)*3 + len(other.bond_types)*3

        data.components = other.data.components

        data.has_substrate = True

        print('hi', data.components)
        print('hello', self.data.components)
        self._film = self.film.copy()
        self._bonds = self.bonds.copy()
        
        already_keys = [i for i in data.components.keys() if isinstance(i, int)]

        offset_to_key = len(already_keys)
        for component_key in self.data.components:
            component = self.data.components.get(component_key)

            if component_key != 'substrate':
                data.components.update({component_key + offset_to_key: component})

        print(other.substrate.shape)
        #old_A = other.substrate_types.get('A')
        #new_A = len(other.types)*2 + len(self.types)*2 + 1
        #old_B = other.substrate_types.get('B')
        #new_B = len(other.types)*2 + len(self.types)*2 + 2
        #print(old_A, new_A, old_B, new_B)
        #other.data.substrate.update({'block_types': {'A': new_A, 'B': new_B}})
        #other.substrate[other.substrate[:,2] == old_A, 2] = new_A
        #other.substrate[other.substrate[:,2] == old_B, 2] = new_B
        print(other.substrate.shape)

        self.film[:,0] += max(other.film[:,0])
        self.film[:,1] += max(other.film[:,1])
        self.film[:,5] += max(other.film[:,5])
        other.substrate[:,0] = np.linspace(1, other.substrate.shape[0], other.substrate.shape[0]).astype(int) + max(self.film[:,0])
        other.substrate[:,1] = max(self.film[:,1]) + 1

        data.atoms = np.row_stack((other.film, self.film, other.substrate))

        self.bonds[:,0] += max(other.bonds[:,0])
        self.bonds[:,2:] += max(other.film[:,0])
        data.bonds = np.row_stack((other.bonds, self.bonds))

        layer = Layer(data)

        # adjust zhi of a resulting film
        data.zhi = np.max(data.atoms[:,5]) + 20

        pairs = []
        pairs = layer._det_pair_film_film(pairs)
        pairs = layer._det_pair_film_substrate(pairs, e_AA=0.745, e_AB=0.745, e_BB=0.745)
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

            if component_key != 'substrate':
                if isinstance(component, list):
                    for sub_component in component:
                        top.append(str(sub_component.get('block_types').get('A')))
                        top.append(str(sub_component.get('block_types').get('B')))

                else:
                    top.append(str(component.get('block_types').get('A')))
                    top.append(str(component.get('block_types').get('B')))

        subs = [str(tp) for tp in layer.substrate_types.values()]
        
        groups.append(f"group\t\t\ttop type {' '.join(top)}")
        groups.append(f"group\t\t\tbottom type {' '.join(bottom)}")
        groups.append(f"group\t\t\tpoly type {' '.join(poly)}")
        groups.append(f"group\t\t\tsubs type {' '.join(subs)}")

        layer.data.templates = t_w_1 + "\n".join(pairs) + "\n".join(groups) + t_w_2 
        layer.data.script_name = "in.welding.txt"
        layer.data.has_substrate = True

        return layer


    def add_substrate(self, substrate_types={'A': 3, 'B': 4}, ratio=0.5):
        """
        add a substrate to the lmp data
        Parameters:
        ----------
        substrate_types (dict): types of beads for the substrate
        
        Returns:
        -------
        N/A
        """

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
        pairs = self._det_pair_film_film(pairs)
        pairs = self._det_pair_film_substrate(pairs)
        pairs.append("\n")

        pairs_mod = []
        pairs_mod = self._det_pair_film_substrate(pairs_mod, e_AA=0.745, e_AB=0.745, e_BB=0.745)
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
        """
        determine types of substrate beads
        """
        rn = np.random.rand(self.num_subsbeads)
        ttype = np.zeros(len(rn))
        
        B = rn > self.substrate_ratio
        A = rn <= self.substrate_ratio
        ttype[B] = self.substrate_types.get('B')
        ttype[A] = self.substrate_types.get('A')
        
        return ttype, sum(A), sum(B)


    def _det_pair_film_film(self, pairs):
        """
        determine all possible interactions between film beads
        Parameters:
        ----------
        pairs (list): list of pairs (normally given as an empty list)
        
        Returns:
        -------
        pairs (list): list of pairs
        """

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
                    type_1 = min(type_Ai, type_Aj)
                    type_2 = max(type_Ai, type_Aj)
                    pairs.append(f"pair_coeff\t\t\t{type_1} {type_2} lj/cut 1.0 1.0 2.5")

                    type_1 = min(type_Ai, type_Bj)
                    type_2 = max(type_Ai, type_Bj)
                    pairs.append(f"pair_coeff\t\t\t{type_1} {type_2} lj/cut 1.0 1.0 2.5")

                    type_1 = min(type_Bi, type_Aj)
                    type_2 = max(type_Bi, type_Aj)
                    pairs.append(f"pair_coeff\t\t\t{type_1} {type_2} lj/cut 1.0 1.0 2.5")

                    type_1 = min(type_Bi, type_Bj)
                    type_2 = max(type_Bi, type_Bj)
                    pairs.append(f"pair_coeff\t\t\t{type_1} {type_2} lj/cut 1.0 1.0 2.5")

        return pairs

    
    def _det_pair_film_substrate(self, pairs, e_AA=None, e_AB=None, e_BB=None):
        """
        determine all possible interactions between film and substrate beads
        Parameters:
        ----------
        pairs (list): list of pairs (normally already filled in with film-film pairs)
        e_AA (float): interaction parameter between A and A
        e_AB (float): interaction parameter between A and B
        e_BB (float): interaction parameter between B and B
        
        Returns:
        -------
        pairs (list): list of pairs
        """

        if e_AA is not None:
            e_AA = e_AA
        else:
            e_AA = 1.0

        if e_AB is not None:
            e_AB = e_AB
        else:
            e_AB = 1.0

        if e_BB is not None:
            e_BB = e_BB
        else:
            e_BB = 1.0

        if e_AB is not None:
            e_AB = e_AB
        else:
            e_AB = 1.0

        for types in self.types:
            type_A = types.get('A')
            type_B = types.get('B')
        
            tp_A = self.substrate_types.get('A')
            tp_B = self.substrate_types.get('B')
            pairs.append(f"pair_coeff\t\t\t{type_A} {tp_A} lj/cut {e_AA} 1.0 2.5")
            pairs.append(f"pair_coeff\t\t\t{type_B} {tp_A} lj/cut {e_AB} 1.0 2.5")
            pairs.append(f"pair_coeff\t\t\t{type_A} {tp_B} lj/cut {e_AB} 1.0 2.5")
            pairs.append(f"pair_coeff\t\t\t{type_B} {tp_B} lj/cut {e_BB} 1.0 2.5")

        return pairs


    def modify_layering(self, from_layering, to_layering):
        """
        modify/fix component_key (layering)
        since the layer ordering may not correctly assigned when reading in an lmp data file
        OR after blended by __init__, component_keys may not be correct

        effectively the same as the Data.modify_component method
        Parameters:
        ----------
        from_c (list): list of current ordering; [0, 1]
        to_c (list): list of new ordering; [1, 0]

        Returns:
        -------
        N/A
        """
        components = self.data.components.copy()

        components_tmp = {}
            
        old_keys = from_layering # ex: [0, 1]
        new_keys = to_layering #ex: [1, 0]

        for old_key, new_key in zip(old_keys, new_keys):
            print(old_key, new_key)
            components_tmp.update({new_key: components.get(old_key)})
            
        components_tmp.update({'substrate': components.get('substrate')})
        print("tmp", components_tmp)
        components = components_tmp

        self.data.components = components


class Script:
    """
    Attributes:
        N/A
    """
    def __init__(self):
        """
        Script constructor
        """
        self.script = []
        self.var = {}
        self.components = {}
        self.substrate = {}
        self.timestep = 0.006


    def set_variables(self, v_type, **kwargs):
        """
        set variables to be used in the lmp script
        Parameters:
        ----------
        v_type (str): either "equal" or "index"
            "equal": cannot be replaced by command line arguments
            "index": can be replaced by command line arguments

        Returns:
        -------
        N/A
        """

        for v_name in kwargs:
            self.var.update({v_name: {"type": v_type, "value": kwargs.get(v_name)}})

        if "alpha" in kwargs:
            self.var.update({
                "epsilon_AA": {"type": "equal", "value": '"1.00 + v_alpha/100"'},
                "epsilon_BB": {"type": "equal", "value": '"1.00 - v_alpha/100"'}
                })

        if "Gamma" in kwargs:
            self.var.update({
                "epsilon_SA": {"type": "equal", "value": '"v_epsilon_AA * v_Gamma + (1 - v_Gamma) * v_epsilon_AB"'},
                "epsilon_SB": {"type": "equal", "value": '"(1 - v_Gamma) * v_epsilon_BB + v_epsilon_AB * v_Gamma"'}
                })


    def set_bc(self, x='p', y='p', z='f'):
        """
        set boundary conditions
        Parameters:
        ----------
        x (str): condition in x axis
        y (str): condition in y axis
        z (str): condition in z axis

        Returns:
        -------
        N/A
        """

        self.bc_x = x
        self.bc_y = y
        self.bc_z = z


    def set_components(self, components):
        """
        set components
        Parameters:
        ----------
        components (dict): components to be used to construct layers and films

        Returns:
        -------
        N/A
        """
        
        self.components.update(components)

        types = []
        bond_types = []
        substrate_types = {}

        for component_key in components:
            component = self.components.get(component_key)

            if component_key != 'substrate':
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

        if 'substrate' in components.keys():
            # substrate

            substrate_types.update(components.get('substrate').get('block_types'))

        self.types =  types
        self.bond_types = bond_types
        self.substrate_types = substrate_types


    def set_timestep(self, timestep):
        """
        set timestep
        Parameters:
        ----------
        timestep (float): timestep

        Returns:
        -------
        N/A
        """
        self.timestep = timestep


    def _add_line(self, line: str):
        """
        add line in the script list
        """
        self.script.append(line)


    def _convert_var(self, v_name: str) -> str:
        """
        helper function that converts given variables to lines that can be directly used in lmp script
        """
        var_temp = self.var.get(v_name)
        v_type = var_temp.get("type")
        v_value = var_temp.get("value")
        
        return f"variable\t\t\t{v_name} {v_type} {v_value}"


    def save(self, lmp_data_name, output):
        """
        save lmp script
        Parameters:
        ----------
        lmp_data_name (str): lmp input data "file name"
        output (str): name of lmp script output

        Returns:
        -------
        N/A
        """

        self._add_line("# LAMMPS script")
        self._add_line("\n## Variables")
        for var in self.var:
            self._add_line(self._convert_var(var))

        self._add_line("\n## Initialization")
        self._add_line("units\t\t\t\tlj")
        self._add_line(f"boundary\t\t\t{self.bc_x} {self.bc_y} {self.bc_z}")
        self._add_line("atom_style\t\t\tbond")

        self._add_line("\n## Neighbor, Forcefield, and Data")
        self._add_line("neighbor\t\t\t0.3 bin")
        self._add_line("neigh_modify\t\tdelay 0 one 2000 page 20000")
        self._add_line("bond_style\t\t\tfene")
        self._add_line("pair_style\t\t\thybrid lj/cut 2.5")
        self._add_line("pair_modify\t\t\tshift yes")

        self._add_line("if \"${sim} == 0\" then &")
        self._add_line(f"  \"read_data {lmp_data_name}\" &")
        self._add_line("  \"#read_restart restart_preequil\" &")
        self._add_line("else &")
        self._add_line("  \"variable index equal 'v_sim - 1'\" &")
        self._add_line("  \"read_restart restart_equil_${index}\" &")
        self._add_line("  \"#read_data data_equil_${index}\" &")


        self._add_line("")
        self._add_line("special_bonds\t\tfene angle no dihedral no lj/coul 0 1 1")
        self._add_line("")
        for ind, types in enumerate(self.bond_types):
            self._add_line(f"bond_coeff\t\t\t{types.get('AA')} 30.0 1.5 ${{epsilon_AA}} 1.0")
            self._add_line(f"bond_coeff\t\t\t{types.get('AB')} 30.0 1.5 ${{epsilon_AB}} 1.0")
            self._add_line(f"bond_coeff\t\t\t{types.get('BB')} 30.0 1.5 ${{epsilon_BB}} 1.0")

        self._add_line("")
        self._add_line("pair_coeff\t\t\t* * none")


        for i, types_i in enumerate(self.types):
            type_Ai = types_i.get('A')
            type_Bi = types_i.get('B')
        
            self._add_line(f"pair_coeff\t\t\t{type_Ai} {type_Ai} lj/cut ${{epsilon_AA}} 1.0 2.5")
            self._add_line(f"pair_coeff\t\t\t{type_Ai} {type_Bi} lj/cut ${{epsilon_AB}} 1.0 2.5")
            self._add_line(f"pair_coeff\t\t\t{type_Bi} {type_Bi} lj/cut ${{epsilon_BB}} 1.0 2.5")

            for j, types_j in enumerate(self.types):
                type_Aj = types_j.get('A')
                type_Bj = types_j.get('B')

                if j > i:
                    type_1 = min(type_Ai, type_Aj)
                    type_2 = max(type_Ai, type_Aj)
                    self._add_line(f"pair_coeff\t\t\t{type_1} {type_2} lj/cut ${{epsilon_AA}} 1.0 2.5")

                    type_1 = min(type_Ai, type_Bj)
                    type_2 = max(type_Ai, type_Bj)
                    self._add_line(f"pair_coeff\t\t\t{type_1} {type_2} lj/cut ${{epsilon_AB}} 1.0 2.5")

                    type_1 = min(type_Bi, type_Aj)
                    type_2 = max(type_Bi, type_Aj)
                    self._add_line(f"pair_coeff\t\t\t{type_1} {type_2} lj/cut ${{epsilon_AB}} 1.0 2.5")

                    type_1 = min(type_Bi, type_Bj)
                    type_2 = max(type_Bi, type_Bj)
                    self._add_line(f"pair_coeff\t\t\t{type_1} {type_2} lj/cut ${{epsilon_BB}} 1.0 2.5")


        self._add_line("")
        for types in self.types:
            type_A = types.get('A')
            type_B = types.get('B')
        
            tp_A = self.substrate_types.get('A')
            tp_B = self.substrate_types.get('B')
            self._add_line(f"pair_coeff\t\t\t{type_A} {tp_A} lj/cut ${{epsilon_SA}} 1.0 2.5")
            self._add_line(f"pair_coeff\t\t\t{type_B} {tp_A} lj/cut ${{epsilon_SB}} 1.0 2.5")
            self._add_line(f"pair_coeff\t\t\t{type_A} {tp_B} lj/cut ${{epsilon_SA}} 1.0 2.5")
            self._add_line(f"pair_coeff\t\t\t{type_B} {tp_B} lj/cut ${{epsilon_SB}} 1.0 2.5")
        
        self._add_line("")
        poly = []
        for types in self.types:
            poly.append(str(types.get('A')))
            poly.append(str(types.get('B')))

        subs = [str(tp) for tp in self.substrate_types.values()]

        self._add_line(f"group\t\t\t\tpoly type {' '.join(poly)}")
        self._add_line(f"group\t\t\t\tsubs type {' '.join(subs)}")

        self._add_line("")
        self._add_line("comm_style\t\t\ttiled")
        self._add_line("#####################################################")
        self._add_line("reset_timestep\t\t0")
        self._add_line(f"if \"${{sim}} == 0\" then &")
        self._add_line(f"  \"velocity    all create ${{T_start}} 900531 dist gaussian\"")

        self._add_line("")
        self._add_line(f"timestep\t\t\t{self.timestep}")

        self._add_line("fix\t\t\t\t\tzwall_lo poly wall/reflect zlo EDGE")
        self._add_line("fix\t\t\t\t\tzwall_hi poly wall/reflect zhi EDGE")
        self._add_line("fix\t\t\t\t\tbal all balance 100000 1.1 rcb out balancing")
        self._add_line("fix\t\t\t\t\tzeroforce subs setforce 0 0 0")
        self._add_line("fix\t\t\t\t\tzeromomentum poly momentum 100 linear 1 1 1")

        self._add_line("")
        self._add_line(f"dump\t\t\t\tsimple all custom 10000 dump.${{sim}}.* id mol type x y z")

        self._add_line(f"fix\t\t\t\t\t1 poly nvt temp ${{T_start}} ${{T_end}} ${{T_damp}}")

        self._add_line("thermo_style\t\tcustom step temp press ke pe epair ebond pxx pyy pzz vol density")
        self._add_line("thermo\t\t\t\t100 ")

        self._add_line("run\t\t\t\t\t10000000")

        self._add_line("unfix\t\t\t\t1")
        self._add_line("unfix\t\t\t\tzeromomentum")

        self._add_line("")
        self._add_line(f"write_data\t\t\tdata_equil_${{sim}}")
        self._add_line(f"write_restart\t\trestart_equil_${{sim}}")

        self._add_line("print \"All done\"")

        script = "\n".join(self.script)
        print(script)

        f = open(f"in.{output}.txt", "w")
        f.write(script)
        f.close()







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
