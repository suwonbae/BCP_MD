import os

"""
TODO: deal with the rest, decorator
"""
def main ():

    #components = {0: {'N': 20, 'f_A': 0.5, 'block_types': {'A': 3, 'B': 4}, 'n': 1800, 'bond_types': {'AA': 4, 'AB': 6, 'BB': 5}}, 'substrate': {'N': 5928, 'f_A': 0.40064102564102566, 'block_types': {'A': 5, 'B': 6}, 'n': 1}, 1: {'N': 20, 'f_A': 0.25, 'block_types': {'A': 1, 'B': 2}, 'n': 1800, 'bond_types': {'AA': 1, 'AB': 3, 'BB': 2}}}

    components = {
            0: [{'N': 20, 'f_A': 0.5, 'block_types': {'A': 3, 'B': 4}, 'n': 1800, 'bond_types': {'AA': 4, 'AB': 6, 'BB': 5}},
                {'N': 20, 'f_A': 0.25, 'block_types': {'A': 1, 'B': 2}, 'n': 1800, 'bond_types': {'AA': 1, 'AB': 3, 'BB': 2}}],
            'substrate': {'N': 5928, 'f_A': 0.40064102564102566, 'block_types': {'A': 5, 'B': 6}, 'n': 1},
            }

    script = LMPscript()
    script.addVariables("index", sim=0, T_start=1.2, T_end=1.2, T_damp=10.0)
    script.addVariables("equal", alpha=1.0, epsilon_AB=0.5, Gamma=0.4)

    script.addComponents(components)
    #script.addSubstrate(types = [3, 4])
    script.set_timestep(0.006)

    script.setBC(x='p', y='p', z='f')

    script.saveScript("test")

class LMPscript:

    def __init__(self):
        self.script = []
        self.var = {}
        self.components = {}
        self.substrate = {}
        self.timestep = 0.006

    def addVariables(self, v_type, **kwargs):

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
            
    def setBC(self, x='p', y='p', z='f'):
        self.bc_x = x
        self.bc_y = y
        self.bc_z = z

    def addComponents(self, components):
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
        self.timestep = timestep

    def _addLine(self, line: str):
        self.script.append(line)

    def _convertVar(self, v_name: str) -> str:

        var_temp = self.var.get(v_name)
        v_type = var_temp.get("type")
        v_value = var_temp.get("value")
        
        return f"variable\t\t\t{v_name} {v_type} {v_value}"

    def saveScript(self, output):
        
        self._addLine("# LAMMPS script")
        self._addLine("\n## Variables")
        for var in self.var:
            self._addLine(self._convertVar(var))

        self._addLine("\n## Initialization")
        self._addLine("units\t\t\t\tlj")
        self._addLine(f"boundary\t\t\t{self.bc_x} {self.bc_y} {self.bc_z}")
        self._addLine("atom_style\t\t\tbond")

        self._addLine("\n## Neighbor, Forcefield, and Data")
        self._addLine("neighbor\t\t\t0.3 bin")
        self._addLine("neigh_modify\t\tdelay 0 one 2000 page 20000")
        self._addLine("bond_style\t\t\tfene")
        self._addLine("pair_style\t\t\thybrid lj/cut 2.5")
        self._addLine("pair_modify\t\t\tshift yes")

        self._addLine("if \"${sim} == 0\" then &")
        self._addLine("  \"read_data data_preequil\" &")
        self._addLine("  \"#read_restart restart_preequil\" &")
        self._addLine("else &")
        self._addLine("  \"variable index equal 'v_sim - 1'\" &")
        self._addLine("  \"read_start restart_equil_${index}\" &")
        self._addLine("  \"#read_data data_equil_${index}\" &")


        self._addLine("")
        self._addLine("special_bonds\t\tfene angle no dihedral no lj/coul 0 1 1")
        self._addLine("")
        for ind, types in enumerate(self.bond_types):
            self._addLine(f"bond_coeff\t\t\t{types.get('AA')} 30.0 1.5 ${{epsilon_AA}} 1.0")
            self._addLine(f"bond_coeff\t\t\t{types.get('AB')} 30.0 1.5 ${{epsilon_AB}} 1.0")
            self._addLine(f"bond_coeff\t\t\t{types.get('BB')} 30.0 1.5 ${{epsilon_BB}} 1.0")

        self._addLine("")
        self._addLine("pair_coeff\t\t\t* * none")


        for i, types_i in enumerate(self.types):
            type_Ai = types_i.get('A')
            type_Bi = types_i.get('B')
        
            self._addLine(f"pair_coeff\t\t\t{type_Ai} {type_Ai} lj/cut ${{epsilon_AA}} 1.0 2.5")
            self._addLine(f"pair_coeff\t\t\t{type_Ai} {type_Bi} lj/cut ${{epsilon_AB}} 1.0 2.5")
            self._addLine(f"pair_coeff\t\t\t{type_Bi} {type_Bi} lj/cut ${{epsilon_BB}} 1.0 2.5")

            for j, types_j in enumerate(self.types):
                type_Aj = types_j.get('A')
                type_Bj = types_j.get('B')

                if j > i:
                    type_1 = min(type_Ai, type_Aj)
                    type_2 = max(type_Ai, type_Aj)
                    self._addLine(f"pair_coeff\t\t\t{type_1} {type_2} lj/cut ${{epsilon_AA}} 1.0 2.5")

                    type_1 = min(type_Ai, type_Bj)
                    type_2 = max(type_Ai, type_Bj)
                    self._addLine(f"pair_coeff\t\t\t{type_1} {type_2} lj/cut ${{epsilon_AB}} 1.0 2.5")

                    type_1 = min(type_Bi, type_Aj)
                    type_2 = max(type_Bi, type_Aj)
                    self._addLine(f"pair_coeff\t\t\t{type_1} {type_2} lj/cut ${{epsilon_AB}} 1.0 2.5")

                    type_1 = min(type_Bi, type_Bj)
                    type_2 = max(type_Bi, type_Bj)
                    self._addLine(f"pair_coeff\t\t\t{type_1} {type_2} lj/cut ${{epsilon_BB}} 1.0 2.5")


        self._addLine("")
        for types in self.types:
            type_A = types.get('A')
            type_B = types.get('B')
        
            tp_A = self.substrate_types.get('A')
            tp_B = self.substrate_types.get('B')
            self._addLine(f"pair_coeff\t\t\t{type_A} {tp_A} lj/cut ${{epsilon_SA}} 1.0 2.5")
            self._addLine(f"pair_coeff\t\t\t{type_B} {tp_A} lj/cut ${{epsilon_SB}} 1.0 2.5")
            self._addLine(f"pair_coeff\t\t\t{type_A} {tp_B} lj/cut ${{epsilon_SA}} 1.0 2.5")
            self._addLine(f"pair_coeff\t\t\t{type_B} {tp_B} lj/cut ${{epsilon_SB}} 1.0 2.5")
        
        self._addLine("")
        poly = []
        for types in self.types:
            poly.append(str(types.get('A')))
            poly.append(str(types.get('B')))

        subs = [str(tp) for tp in self.substrate_types.values()]

        self._addLine(f"group\t\t\t\tpoly type {' '.join(poly)}")
        self._addLine(f"group\t\t\t\tsubs type {' '.join(subs)}")

        self._addLine("")
        self._addLine("comm_style\t\t\ttiled")
        self._addLine("#####################################################")
        self._addLine("reset_timestep\t\t0")
        self._addLine(f"if \"${{sim}} == 0\" then &")
        self._addLine(f"  \"velocity    all create ${{T_start}} 900531 dist gaussian\"")

        self._addLine("")
        self._addLine(f"timestep\t\t\t{self.timestep}")

        self._addLine("fix\t\t\t\t\tzwall_lo poly wall/reflect zlo EDGE")
        self._addLine("fix\t\t\t\t\tzwall_hi poly wall/reflect zhi EDGE")
        self._addLine("fix\t\t\t\t\tbal all balance 100000 1.1 rcb out balancing")
        self._addLine("fix\t\t\t\t\tzeroforce subs setforce 0 0 0")
        self._addLine("fix\t\t\t\t\tzeromomentum poly momentum 100 linear 1 1 1")

        self._addLine("")
        self._addLine(f"dump\t\t\t\tsimple all custom 10000 dump.${{sim}}.* id mol type x y z")

        self._addLine(f"fix\t\t\t\t\t1 poly nvt temp ${{T_start}} ${{T_end}} ${{T_damp}}")

        self._addLine("thermo_style\t\tcustom step temp press ke pe epair ebond pxx pyy pzz vol density")
        self._addLine("thermo\t\t\t\t100 ")

        self._addLine("run\t\t\t\t\t10000000")

        self._addLine("unfix\t\t\t\t1")
        self._addLine("unfix\t\t\t\tzeromomentum")

        self._addLine("")
        self._addLine(f"write_data\t\t\tdata_equil_${{sim}}")
        self._addLine(f"write_restart\t\trestart_equil_${{sim}}")

        self._addLine("print \"All done\"")

        script = "\n".join(self.script)
        print(script)

        f = open(f"in.{output}.txt", "w")
        f.write(script)
        f.close()

if __name__ == "__main__":
    main()
       
