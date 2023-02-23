import sqlite3
import numpy as np
import re
import os
import matplotlib.pyplot as plt

__all__ = ['Sqldb', 'parse', 'parse_log', 'parse_path']

class Sqldb:

    def __init__(self, db):

        try:
            conn = sqlite3.connect(db)

            print("Connection is established")

            self.conn = conn

            self.activate_foreign_keys()

        except sqlite3.Error:

            print(sqlite3.Error)


    def activate_foreign_keys(self):
        
        conn = self.conn
        c = conn.cursor()
        c.execute("PRAGMA foreign_keys = ON;")
        print("## foreign_keys = ON\n")


    def deactivate_foriegn_keys(self):

        conn = self.conn
        c = conn.cursor()
        c.execute("PRAGMA foreign_keys = OFF;")
        print("## foreign_keys = OFF\n")        


    def identify_chain(self, components):

        conn = self.conn
        c = conn.cursor()

        query = """
        SELECT chain_id
        FROM chains
        WHERE N = ? and f_A = ?;
        """

        print("## chain")
        for component_key in components:
            if not isinstance(component_key, int): continue

            component = components.get(component_key)

            if isinstance(component, list):

                for sub_component in component:

                    N = int(sub_component.get('N'))
                    f_A = float(sub_component.get('f_A'))

                    tbl = c.execute(query, (N, f_A)).fetchall()
                    if tbl:
                        sub_component.update({"chain_id": tbl[0][0]})
                        print(f"# found chain_id: {tbl[0][0]}")
                    else:
                        c.execute("INSERT INTO chains (N, f_A) VALUES (?, ?);", (N, f_A))
                        sub_component.update({"chain_id": c.lastrowid})
                        print(f"# added a new chain (chain_id: {c.lastrowid})")

            elif isinstance(component, dict):

                N = int(component.get('N'))
                f_A = float(component.get('f_A'))

                tbl = c.execute(query, (N, f_A)).fetchall()
                if tbl:
                    component.update({"chain_id": tbl[0][0]})
                    print(f"# found chain_id: {tbl[0][0]}")
                else:
                    c.execute("INSERT INTO chains (N, f_A) VALUES (?, ?);", (N, f_A))
                    component.update({"chain_id": c.lastrowid})
                    print(f"# added a new chain (chain_id: {c.lastrowid})")
            
            components.update({component_key: component})

        print("")
        conn.commit()

        return components


    def identify_layer(self, components):

        conn = self.conn
        c = conn.cursor()

        print("## layer")
        for component_key in components:
            if not isinstance(component_key, int): continue

            component = components.get(component_key)

            bc_x, bc_y, bc_z = 'p', 'p', 'f'

            if isinstance(component, list):

                for sub_component in component:
                    if "bc_x" in sub_component:
                        bc_x = component.get('bc_x')
                    if "bc_y" in sub_component:
                        bc_y = component.get('bc_y')
                    if "bc_z" in sub_component:
                        bc_z = component.get('bc_z')

                    if 'presim_id' in sub_component:
                        presim = True

                        query = """
                        SELECT layer_id
                        FROM layers
                        WHERE chain_id = ? AND n = ? AND presim_id = ? AND bc_x = ? AND bc_y = ? AND bc_z = ?;
                        """

                        tbl = c.execute(query, (sub_component.get('chain_id'), sub_component.get('n'), sub_component.get('presim_id'), bc_x, bc_y, bc_z)).fetchall()
                    else:
                        presim = False

                        query = """
                        SELECT layer_id
                        FROM layers
                        WHERE chain_id = ? AND n = ? AND bc_x = ? AND bc_y = ? AND bc_z = ?;
                        """                

                        tbl = c.execute(query, (sub_component.get('chain_id'), sub_component.get('n'), bc_x, bc_y, bc_z)).fetchall()

                    if tbl:
                        sub_component.update({"layer_id": tbl[0][0]})
                        print(f"# found layer_id: {tbl[0][0]}")
                    else:
                        if presim:
                            c.execute("INSERT INTO layers (chain_id, n, presim_id, bc_x, bc_y, bc_z) VALUES (?, ?, ?, ?, ?, ?);", (sub_component.get('chain_id'), sub_component.get('n'), sub_component.get('presim_id'), bc_x, bc_y, bc_z))
                            sub_component.update({"layer_id": c.lastrowid})
                        else:
                            c.execute("INSERT INTO layers (chain_id, n, bc_x, bc_y, bc_z) VALUES (?, ?, ?, ?, ?);", (sub_component.get('chain_id'), sub_component.get('n'), bc_x, bc_y, bc_z))
                            sub_component.update({"layer_id": c.lastrowid})
                        print(f"# added a new layer (layer_id: {c.lastrowid})")

            elif isinstance(component, dict):

                if "bc_x" in component:
                    bc_x = component.get('bc_x')
                if "bc_y" in component:
                    bc_y = component.get('bc_y')
                if "bc_z" in component:
                    bc_z = component.get('bc_z')

                if 'presim_id' in component:
                    presim = True

                    query = """
                    SELECT layer_id
                    FROM layers
                    WHERE chain_id = ? AND n = ? AND presim_id = ? AND bc_x = ? AND bc_y = ? AND bc_z = ?;
                    """

                    tbl = c.execute(query, (component.get('chain_id'), component.get('n'), component.get('presim_id'), bc_x, bc_y, bc_z)).fetchall()
                else:
                    presim = False

                    query = """
                    SELECT layer_id
                    FROM layers
                    WHERE chain_id = ? AND n = ? AND bc_x = ? AND bc_y = ? AND bc_z = ?;
                    """                

                    tbl = c.execute(query, (component.get('chain_id'), component.get('n'), bc_x, bc_y, bc_z)).fetchall()

                if tbl:
                    component.update({"layer_id": tbl[0][0]})
                    print(f"# found layer_id: {tbl[0][0]}")
                else:
                    if presim:
                        c.execute("INSERT INTO layers (chain_id, n, presim_id, bc_x, bc_y, bc_z) VALUES (?, ?, ?, ?, ?, ?);", (component.get('chain_id'), component.get('n'), component.get('presim_id'), bc_x, bc_y, bc_z))
                        component.update({"layer_id": c.lastrowid})
                    else:
                        c.execute("INSERT INTO layers (chain_id, n, bc_x, bc_y, bc_z) VALUES (?, ?, ?, ?, ?);", (component.get('chain_id'), component.get('n'), bc_x, bc_y, bc_z))
                        component.update({"layer_id": c.lastrowid})
                    print(f"# added a new layer (layer_id: {c.lastrowid})")

            components.update({component_key: component})

        print("")
        conn.commit()

        return components


    def identify_film(self, components, layering=None, comment=None):

        conn = self.conn
        c = conn.cursor()

        # if component is list: blend, ditc: pure
        # 1 component -> single layer
        # more than 1 -> multi-layered

        if layering is not None:
            components_tmp = {}
            
            old_keys = list(components.keys()) # ex: [0, 1]
            new_keys = layering #ex: [1, 0]

            for old_key, new_key in zip(old_keys, new_keys):
                print(old_key, new_key)
                components_tmp.update({new_key: components.get(old_key)})
            
            print("tmp", components_tmp)
            components = components_tmp

        query = """
        SELECT film_id
        FROM film_info
        WHERE layer_id = ? AND layering = ? AND A = ? AND B = ? AND AA = ? AND AB = ? AND BB= ?;
        """
        
        if len(components) > 1:
            print("multi-layered")

            tbls = []
            for component_key in components:
                if not isinstance(component_key, int): continue

                component = components.get(component_key)

                if isinstance(component, list):
                    print("blend")
                    
                    for sub_component in component:
                        type_A = int(sub_component.get('block_types').get('A'))
                        type_B = int(sub_component.get('block_types').get('B'))
                        type_AA = int(sub_component.get('bond_types').get('AA'))
                        type_AB = int(sub_component.get('bond_types').get('AB'))
                        type_BB = int(sub_component.get('bond_types').get('BB'))

                        tbl = c.execute(query, (sub_component.get('layer_id'), component_key, type_A, type_B, type_AA, type_AB, type_BB)).fetchall()
                        tbls.append(tbl)

                elif isinstance(component, dict):
                    print("pure")

                    type_A = int(component.get('block_types').get('A'))
                    type_B = int(component.get('block_types').get('B'))
                    type_AA = int(component.get('bond_types').get('AA'))
                    type_AB = int(component.get('bond_types').get('AB'))
                    type_BB = int(component.get('bond_types').get('BB'))

                    tbl = c.execute(query, (component.get('layer_id'), component_key, type_A, type_B, type_AA, type_AB, type_BB)).fetchall()
                    tbls.append(tbl)

            for ind in range(len(tbls)):
                if ind == 0:
                    tbl = tbls[ind]
                else:
                    tbl = list(set(tbl) & set(tbls[ind]))

            if tbl:
                film_id = tbl[0][0]
                components.update({'film_id': film_id})
                print(f"# found film_id: {film_id}")
            else:
                c.execute("INSERT INTO films (comment) VALUES (?);", (comment, ))
                film_id = c.lastrowid

                for component_key in components:
                    if not isinstance(component_key, int): continue

                    component = components.get(component_key)

                    if isinstance(component, list):
                        print("blend")
                        for sub_component in component:
                        
                            type_A = int(sub_component.get('block_types').get('A'))
                            type_B = int(sub_component.get('block_types').get('B'))
                            type_AA = int(sub_component.get('bond_types').get('AA'))
                            type_AB = int(sub_component.get('bond_types').get('AB'))
                            type_BB = int(sub_component.get('bond_types').get('BB'))
                        
                            c.execute("INSERT INTO film_info (film_id, layer_id, layering, A, B, AA, AB, BB) VALUES (?, ?, ?, ?, ?, ?, ?, ?);", (film_id, sub_component.get('layer_id'), component_key, type_A, type_B, type_AA, type_AB, type_BB))

                    elif isinstance(component, dict):
                        print("pure")

                        type_A = int(component.get('block_types').get('A'))
                        type_B = int(component.get('block_types').get('B'))
                        type_AA = int(component.get('bond_types').get('AA'))
                        type_AB = int(component.get('bond_types').get('AB'))
                        type_BB = int(component.get('bond_types').get('BB'))
                        
                        c.execute("INSERT INTO film_info (film_id, layer_id, layering, A, B, AA, AB, BB) VALUES (?, ?, ?, ?, ?, ?, ?, ?);", (film_id, component.get('layer_id'), component_key, type_A, type_B, type_AA, type_AB, type_BB))
                
                components.update({'film_id': film_id})
                print(f"# added a new film (film_id: {film_id})")

        elif len(components) == 1:
            print("single layer")

            component_key = list(components.keys())[0]
            component = components.get(component_key)

            if isinstance(component, list):
                print("blend")
                
                tbls = []
                for sub_component in component:
                    type_A = int(sub_component.get('block_types').get('A'))
                    type_B = int(sub_component.get('block_types').get('B'))
                    type_AA = int(sub_component.get('bond_types').get('AA'))
                    type_AB = int(sub_component.get('bond_types').get('AB'))
                    type_BB = int(sub_component.get('bond_types').get('BB'))

                    tbl = c.execute(query, (sub_component.get('layer_id'), component_key, type_A, type_B, type_AA, type_AB, type_BB)).fetchall()
                    tbls.append(tbl)

                for ind in range(len(component)):
                    if ind == 0:
                        tbl = tbls[ind]
                    else:
                        tbl = list(set(tbl) & set(tbls[ind]))

                if tbl:
                    film_id = tbl[0][0]
                    component.append({'film_id': film_id})
                    print(f"# found film_id: {film_id}")
                else:
                    c.execute("INSERT INTO films (comment) VALUES (?);", (comment, ))
                    film_id = c.lastrowid
                    for sub_component in component:
                        type_A = int(sub_component.get('block_types').get('A'))
                        type_B = int(sub_component.get('block_types').get('B'))
                        type_AA = int(sub_component.get('bond_types').get('AA'))
                        type_AB = int(sub_component.get('bond_types').get('AB'))
                        type_BB = int(sub_component.get('bond_types').get('BB'))

                        c.execute("INSERT INTO film_info (film_id, layer_id, layering, A, B, AA, AB, BB) VALUES (?, ?, ?, ?, ?, ?, ?, ?);", (film_id, sub_component.get('layer_id'), component_key, type_A, type_B, type_AA, type_AB, type_BB))
                    
                    component.append({'film_id': film_id})
                    print(f"# added a new film (film_id: {film_id})")

            elif isinstance(component, dict):
                print("pure")
                
                type_A = int(component.get('block_types').get('A'))
                type_B = int(component.get('block_types').get('B'))
                type_AA = int(component.get('bond_types').get('AA'))
                type_AB = int(component.get('bond_types').get('AB'))
                type_BB = int(component.get('bond_types').get('BB'))

                tbl = c.execute(query, (component.get('layer_id'), component_key, type_A, type_B, type_AA, type_AB, type_BB)).fetchall()

                if tbl:
                    film_id = tbl[0][0]
                    component.update({"film_id": film_id})
                    print(f"# found film_id: {film_id}")
                else:
                    c.execute("INSERT INTO films (comment) VALUES (?);", (comment, ))
                    film_id = c.lastrowid
                    c.execute("INSERT INTO film_info (film_id, layer_id, layering, A, B, AA, AB, BB) VALUES (?, ?, ?, ?, ?, ?, ?, ?);", (film_id, component.get('layer_id'), component_key, type_A, type_B, type_AA, type_AB, type_BB))
                    component.update({"film_id": film_id})
                    print(f"# added a new film (film_id: {film_id})")                    

            components.update({component_key: component})
        
        print("")
        conn.commit()
    
        return components


    def add_sim(self, components, path, substrate_id=None, comment=None):

        conn = self.conn
        c = conn.cursor()

        if substrate_id is None:
            c.execute("INSERT INTO sims (film_id, comment, path) VALUES (?, ?, ?);", (components.get('film_id'), comment, path))
        else:
            c.execute("INSERT INTO sims (film_id, substrate_id, comment, path) VALUES (?, ?, ?, ?);", (components.get('film_id'), substrate_id, comment, path))
        sim_id = c.lastrowid
        self.sim_id = sim_id

        conn.commit()


    def add_sim_dynamics(self, sequence, ensemble='nvt', T_start=1.2, T_end=1.2, T_damp=10.0, timestep=0.006, steps=10000000, log=None, comment=None):
        conn = self.conn
        c = conn.cursor()

        c.execute("""
        INSERT INTO sim_dynamics
        (sim_id, sequence, ensemble, T_start, T_end, T_damp, timestep, steps, log, comment)
        VALUES
        (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
        """, (self.sim_id, sequence, ensemble, T_start, T_end, T_damp, timestep, steps, log, comment))

        conn.commit()


    def add_sim_parameters(self, e_AB, alpha=None, Gamma=None, e_AA=None, e_BB=None, e_SA=None, e_SB=None, **kwargs):

        conn = self.conn
        c = conn.cursor()

        c.execute("""
        INSERT INTO sim_parameters (sim_id, e_AB, alpha, Gamma, e_AA, e_BB, e_SA, e_SB)
        VALUES
        (?, ?, ?, ?, ?, ?, ?, ?);""", (self.sim_id, e_AB, alpha, Gamma, e_AA, e_BB, e_SA, e_SB))

        conn.commit()


    def _generate_condition(self, keys, args):

        conditions = []
        values = []

        for key in keys:
            if key in args:
                condition = f"{key} = ?"
                conditions.append(condition)
                values.append(args.get(key))
        
        return conditions, values


    def query_chain(self, **kwargs):

        query = """
        SELECT
            *
        FROM
            chains
        """
        conditions, values = self._generate_condition(['N', 'f_A'], kwargs)

        if len(conditions) > 0:
            query += f"WHERE {' AND '.join(conditions)}\n"
        query += f"ORDER BY N;"
 
        conn = self.conn
        c = conn.cursor()

        if len(conditions) > 0:
            tbl = c.execute(query, tuple(values)).fetchall()
        else:
            tbl = c.execute(query).fetchall()


        return tbl


    def query_layer(self, **kwargs):

        query = """
        SELECT
            *
        FROM
            layers
        """

        conditions, values = self._generate_condition(['layer_id', 'chain_id', 'n'], kwargs)

        if len(conditions) > 0:
            query += f"WHERE {' AND '.join(conditions)}\n"
        query += f"ORDER BY n;"

        conn = self.conn
        c = conn.cursor()

        if len(conditions) > 0:
            tbl = c.execute(query, tuple(values)).fetchall()
        else:
            tbl = c.execute(query).fetchall()


        return tbl


    def query_film_by_layerid(self, layer_id):

        query = """
        SELECT
            film_id,
            CASE
                WHEN pure = 1 AND blend = 1
                    THEN 'single/pure'
                WHEN pure = 2 AND blend = 1
                    THEN 'single/blend'
                WHEN pure = 2 AND blend = 2
                    THEN 'bilayer'
                ELSE 'etc'
            END category,
            N_total,
            used_layer AS layer,
            used_N AS chain_length,
            used_f_A AS f_A,
            used_nn AS n,
            used_presim AS presim
        FROM
            (SELECT 
                film_id,
                count(film_id) AS pure,
                count(DISTINCT(layering)) AS blend,
                sum(layers.n*chains.N) AS N_total,
                group_concat(film_info.layer_id) AS used_layer,
                group_concat(chains.N) AS used_N,
                group_concat(chains.f_A) AS used_f_A,
                group_concat(layers.n) AS used_nn,
                group_concat(presims.comment) AS used_presim
            FROM
                film_info
            JOIN layers
            ON layers.layer_id = film_info.layer_id
            JOIN chains
            ON chains.chain_id = layers.chain_id
            JOIN presims
            ON presims.presim_id = layers.presim_id
            WHERE film_id IN (
                    SELECT film_id
                    FROM film_info
                    WHERE layer_id = ?
                    )
            GROUP BY film_id
            );
        """

        conn = self.conn
        c = conn.cursor()

        tbl = c.execute(query, (layer_id,)).fetchall()

        return tbl


    def query_film_by_filmid(self, film_id):

        query = """
        SELECT
            film_id,
            CASE
                WHEN pure = 1 AND blend = 1
                    THEN 'single/pure'
                WHEN pure = 2 AND blend = 1
                    THEN 'single/blend'
                WHEN pure = 2 AND blend = 2
                    THEN 'bilayer'
                ELSE 'etc'
            END category,
            N_total,
            used_layer AS layer,
            used_N AS chain_length,
            used_f_A AS f_A,
            used_nn AS n,
            used_presim AS presim
        FROM
            (SELECT 
                film_id,
                count(film_id) AS pure,
                count(DISTINCT(layering)) AS blend,
                sum(layers.n*chains.N) AS N_total,
                group_concat(film_info.layer_id) AS used_layer,
                group_concat(chains.N) AS used_N,
                group_concat(chains.f_A) AS used_f_A,
                group_concat(layers.n) AS used_nn,
                group_concat(presims.comment) AS used_presim
            FROM
                film_info
            JOIN layers
            ON layers.layer_id = film_info.layer_id
            JOIN chains
            ON chains.chain_id = layers.chain_id
            JOIN presims
            ON presims.presim_id = layers.presim_id
            WHERE film_id = ?
            GROUP BY film_id
            );
        """

        conn = self.conn
        c = conn.cursor()

        tbl = c.execute(query, (film_id,)).fetchall()

        return tbl


    def query_sim(self, **kwargs):

        query = """
        SELECT
            *
        FROM
            sims
        """        

        conditions, values = self._generate_condition(['sim_id', 'film_id'], kwargs)

        if len(conditions) > 0:
            query += f"WHERE {' AND '.join(conditions)}\n"

        conn = self.conn
        c = conn.cursor()

        if len(conditions) > 0:
            tbl = c.execute(query, tuple(values)).fetchall()
        else:
            tbl = c.execute(query).fetchall()


        return tbl


    def report_sim(self, sim_id):

        conn = self.conn
        c = conn.cursor()

        query_simdyn = """
        SELECT
            sequence,
            T_start,
            T_end,
            timestep,
            steps
        FROM
            sim_dynamics
        WHERE sim_id = ?;
        """

        query_simpar = """
        SELECT
            e_AB,
            alpha,
            Gamma
        FROM
            sim_parameters
        WHERE sim_id = ?;
        """

        dyn = c.execute(query_simdyn, (sim_id,)).fetchall()
        par = c.execute(query_simpar, (sim_id,)).fetchall()

        dyn = np.asarray(dyn)
        par = np.asarray(par)
        
        T = []
        t = []
        for dyn_line in dyn:
            sequence = dyn_line[0]
            T_start = dyn_line[1]
            T_end = dyn_line[2]
            timestep = dyn_line[3]
            steps = dyn_line[4]

            T.append(T_start)
            T.append(T_end)
            t.append(sequence*steps*timestep)
            t.append((sequence+1)*steps*timestep)          

        e_AB, alpha, Gamma = par[0]

        print(f"* e_AB = {e_AB}, alpha = {alpha}, Gamma = {Gamma}")
        print(f"* substrate: {self.get_substrate_info(sim_id)}")
        print(f"* t: {t[0]} to {t[-1]}")

        path = self.get_path(sim_id)
        dirs = next(os.walk(path))[1]
        dirs = natural_sort(dirs)
        print(f"* {dirs[0]} to {dirs[-1]}")

        try:
            plt.rc('font', size=13)
            plt.rcParams["text.usetex"] = True
        except:
            print("no tex")

        fig, ax = plt.subplots(figsize=(3,2), dpi=200)
        ax.plot(t, T)
        x_min = min(t)
        x_max = max(t)
        y_min = min(T)
        y_max = max(T)

        offset = 0.1
        xticks = np.linspace(t[0], t[-1], 4)
        xticklabels = ["${{{t}}}$".format(t=t[0]), '', '', "${{{t}}}$".format(t=t[-1])]
        ax.text(x_min, y_max + offset, f"$\\epsilon_{{\mathrm{{AB}}}}={e_AB}$, $\\alpha={alpha}$, $\\Gamma={Gamma}$", ha='left', va='top', fontsize=9)
        ax.set_xlabel('time')
        ax.set_ylabel('temperature')
        ax.set_xlim(x_min, x_max)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        ax.set_ylim(y_min - offset, y_max + offset)
        plt.tight_layout(pad=0.5, h_pad=None, w_pad=None, rect=None)
        plt.savefig(f"sim_id_{sim_id}_dynamics.png", dpi=200)


    def get_types(self, sim_id):

        query = """
        SELECT
            film_info.*
        FROM
            film_info
        JOIN sims
        ON sims.film_id = film_info.film_id
        WHERE sim_id = ?
        ORDER BY layering;
        """

        conn = self.conn
        c = conn.cursor()

        tbl = c.execute(query, (sim_id,)).fetchall()
        res = np.asarray(tbl)

        types = {}

        types.update({'all': list(res[:, 3:5].ravel())})
        for i in range(res.shape[0]):
            types.update({f"C{i}": list(res[i, 3:5])})

        A = res[:,3]
        B = res[:,4]
        AA = res[:,5]
        AB = res[:,6]
        BB = res[:,7]

        types.update({'A': list(A)})
        types.update({'B': list(B)})
        types.update({'AA': list(AA)})
        types.update({'AB': list(AB)})
        types.update({'BB': list(BB)})

        return types


    def get_path(self, sim_id):

        query = """
        SELECT
            path
        FROM
            sims
        WHERE sim_id = ?;
        """

        conn = self.conn
        c = conn.cursor()

        tbl = c.execute(query, (sim_id,)).fetchone()
        path = tbl[0]

        return path


    def get_substrate_info(self, sim_id):
        
        query = """
        SELECT
            substrates.comment
        FROM
            substrates
        JOIN sims
        ON sims.substrate_id = substrates.substrate_id
        WHERE sim_id = ?;
        """

        conn = self.conn
        c = conn.cursor()

        tbl = c.execute(query, (sim_id,)).fetchone()
        comment = tbl[0]

        return comment


    def close(self):
        '''
        close connection
        '''

        conn = self.conn
        conn.close()
        print("Connection is closed")


def parse(path=None, pattern_suffix=None, dynamics={}, parameters={}):
    '''
    call parse_log and parse_path
    '''
    if path is None:
        path = os.getcwd()

    dynamics = parse_log(path, pattern_suffix, dynamics)
    parameters = parse_path(path, parameters)

    for dynamics_key in dynamics:
        
        dynamics_tmp = dynamics.get(dynamics_key)

        if 'alpha' in dynamics_tmp:
            parameters.update({'alpha': dynamics_tmp.get('alpha')})
        
        if 'e_AB' in dynamics_tmp:
            parameters.update({'e_AB': dynamics_tmp.get('e_AB')})

        if 'Gamma' in dynamics_tmp:
            parameters.update({'Gamma': dynamics_tmp.get('Gamma')})

    return dynamics, parameters


def parse_log(path, pattern_suffix=None, dynamics={}):
    '''
    extract from log files located in equil directories in the given path 
    '''

    dir_pattern = re.compile(r"(?<=equil_)\d+$")
    if pattern_suffix is not None:
        dir_pattern = re.compile(r"(?<=equil_)\d+"+f"({pattern_suffix})")

    dirs = os.listdir(path)
    numbers = []

    for directory in dirs:
        match = dir_pattern.search(directory)
        if match:
            if pattern_suffix is not None:
                numbers.append(int(match.group(0).split('_')[0]))
            else:
                numbers.append(int(match.group(0)))

    numbers.sort()

    pattern_alpha = re.compile(r"(?<=alpha equal )[+-]?\d+(?:\.\d+)?")
    pattern_epsilon_AB = re.compile(r"(?<=epsilon_AB equal )[+-]?\d+(?:\.\d+)?")
    pattern_Gamma = re.compile(r"(?<=Gamma equal )[+-]?\d+(?:\.\d+)?")

    pattern_timestep = re.compile("^timestep")
    pattern_steps = re.compile("^run")
    pattern_ensemble = re.compile("nvt temp [+-]?\d+(?:\.\d+)? [+-]?\d+(?:\.\d+)? [+-]?\d+(?:\.\d+)?")

    for number in numbers:
        res = {}

        steps = 0
        directory = f"equil_{number}"
        if pattern_suffix is not None:
            directory += pattern_suffix
        dir_path = os.path.join(path, directory)
        if not os.path.exists(os.path.join(dir_path, "log.lammps")):
            print(f"{directory} is missing log.lammps")
        else:
            f = open(os.path.join(dir_path, "log.lammps"), "r")
            flag = 0
            count = 0
            while (flag == 0):
                line = f.readline()
                count += 1
                if len(line) > 0:

                    match = pattern_alpha.search(line)
                    if match:
                        res.update({'alpha': float(match.group(0))})

                    match = pattern_epsilon_AB.search(line)
                    if match:
                        res.update({'e_AB': float(match.group(0))})

                    match = pattern_Gamma.search(line)
                    if match:
                        res.update({'Gamma': float(match.group(0))})

                    match = pattern_timestep.search(line)
                    if match:
                        res.update({"timestep": float(line.split()[1])})

                    match = pattern_steps.search(line)
                    if match:
                        print(f"{number} {int(line.split()[1])}")
                        steps += int(line.split()[1])
                        #res.update({"steps": int(line.split()[1])})

                    match = pattern_ensemble.search(line)
                    if match:
                        segs = match.group(0).split()
                        res.update({"ensemble": segs[0],
                            "T_start": float(segs[2]),
                            "T_end": float(segs[3]),
                            "T_damp": float(segs[4])
                            })

                    #if count > 250:
                    #    flag = 1
                else:
                    flag = 1

            f.close()
            res.update({'steps': steps})
            
        dynamics.update({number: res})

    return dynamics


def parse_path(path, parameters={}):
    '''
    extract parameters from the given path
    '''

    pattern_eAB = re.compile(r"(?<=eAB_)[+-]?\d+(?:\.\d+)?")
    match = pattern_eAB.search(path)
    if match:
        parameters.update({"e_AB": float(match.group(0))})

    pattern_alpha = re.compile(r"(?<=alpha_)[+-]?\d+(?:\.\d+)?")
    match = pattern_alpha.search(path)
    if match:
        parameters.update({"alpha": float(match.group(0))})

    pattern_T = re.compile(r"(?<=T_)[+-]?\d+(?:\.\d+)?")
    match = pattern_T.search(path)
    if match:
        parameters.update({"T": float(match.group(0))})

    pattern_Gamma = re.compile(r"(?<=Gamma_)[+-]?\d+(?:\.\d+)?")
    match = pattern_Gamma.search(path)
    if match:
        parameters.update({"Gamma": float(match.group(0))})

    return parameters


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]

    return sorted(l, key=alphanum_key)