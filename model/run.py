import os
from genLayeredFilm import *

def main():

    if not os.path.exists('results'):
        os.makedirs('results')

    '''
    If the corresponding data file exists, read it in.
    Otherwise, call a randomwalk function to generate a data file.

    In this, I will specify what types of beads construct polymer chains and substrates,
    meaning that I do not have to keep mentioning bead types every time.

    However, when I read in a data file, it is not easy to tell right away identify whether
    the model has a substrate or not (Maybe can have this info in every data file or in DB).

    When generating a randomwalk model, I need the following information:
    how many components, N, f_A, A, B, n, xlo, xhi, ylo, yhi, zlo, zhi, bc_x, bc_y, bc_z
    
    model = {
        components: {
            0: {"N": 20, "f_A": 0.25, "n": 4800, "type_A": 1, "type_B": 2},
            #1: {"N": 22, "f_A": 0.5, "n": 4364, "type_A": 3, "type_B": 4},
        },
        substrate: {},
        box: {"xlo": 0.0, "xhi": 90.0, "ylo": 0.0, "yhi": 57.0, "zlo": 0, "zhi": 20.0},
        bc: {"x": "p", "y": "p", "z": "f"}
    }
    '''

    '''
    1. pure bcp film
    '''

    components = {
            0: {"N": 20, "f_A": 0.25, "block_types": {'A': 1, 'B': 2}, "bond_types": {'AA': 1, 'BB': 2, 'AB': 3}, "n": 2400},
            }

    # RandomWalk(components) -> data file
    # read in data (by self-avoiding random walk)
    #data = LammpsData("data_C20n2400_wo_substate.txt")
    #data.probe()

    # generate a layer
    #layer = Layer(data)

    # add a substrate
    #layer.add_substrate(substrate_types={'A': 3, 'B':4}, ratio=0.5)
    
    # write the new data
    #layer.data.write_data("data_tuned.txt", header="lmp data file")

    '''
    2. bilayer film
    '''
    # read in film(s) (already with their own substrates)
    #data_C = LammpsData("data_C20n4800.txt")
    #data_C.probe()
    #layer_top = Layer(data_C)
    #layer_bottom = Layer(data_C, {0: {'block_types': {'A': 3, 'B': 4}, 'bond_types': {'AA': 4, 'BB': 5, 'AB': 6}}})

    #data_C = LammpsData("data_C20n4800.txt")
    #data_CC = LammpsData("data_tuned.txt")
    #data_C.probe()
    #data_CC.probe()
    #layer_top = Layer(data_C)
    #layer_bottom = Layer(data_CC)

    #layer = layer_top / layer_bottom
    
    #layer.data.write_data("test.txt", header="lmp data file")
    
    '''
    pure
    '''
    #data = LammpsData("data_pure_1.txt")
    #print(data.components)
    #layer = Layer(data)
    #print(f"before adding substrate: has substrate {layer.data.has_substrate}")
    #if not layer.data.has_substrate:
    #    layer.add_substrate(substrate_types={'A': 3, 'B': 4})
    #print(f"after adding substrate: has substrate {layer.data.has_substrate}")
    #layer.data.write("data_pure_1_result.txt")

    import time
    '''
    pure bilayer
    '''
    #data_C = LammpsData("data_pure_3.txt")
    #data_L = LammpsData("data_pure_2.txt")
    #print(data_C.components)
    #print(data_L.components)

    #layer_C = Layer(data_C)
    #layer_L = Layer(data_L)

    #new_layer = layer_C / layer_L
    #new_layer.data.probe()
    #new_layer.data.tailor(0, block_types={'A': 3, 'B': 4}, bond_types={'AA': 4, 'AB': 6, 'BB': 5})
    #new_layer.data.tailor(1, block_types={'A': 1, 'B': 2}, bond_types={'AA': 1, 'AB': 3, 'BB': 2})
    #print(new_layer.data.components)

    '''
    blend
    '''
    # read in a data file of a blend
    #data = LammpsData("data_blend_1.txt")
    #print("before tailoring")
    #print(data.components)
    # fix bond_types of component 1
    #data.tailor(1, bond_types={'AA': 4, 'AB': 6, 'BB': 5})
    #print("after tailroing")
    #print(data.components)
    #data.write("data_blend_1_result.txt")


    '''
    a blend on top of a blend
    '''
    data_1 = LammpsData("data_blend_2.txt")
    data_2 = LammpsData("data_blend_2.txt")
    print("before layering")
    print(data_1.components)
    print(data_2.components)
    layer_1 = Layer(data_1, blend=[0, 1])
    layer_2 = Layer(data_2, blend=[0, 1])
    new_layer = layer_1 / layer_2
    print("after layering")
    print(new_layer.data.components)
    # for blend layer.data.tailor does not work
    #new_layer.data.tailor(0, block_types={'A': 3, 'B': 4}, bond_types={'AA': 4, 'AB': 6, 'BB': 5})
    new_layer.data.write("data_blend_2_result.txt")

if __name__ == '__main__':
    main()
