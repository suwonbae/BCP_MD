import numpy as np
from matplotlib.colors import ListedColormap

cmaps = {}
    
W2B = np.dstack((np.linspace(1,0,256), np.linspace(1,0,256), np.linspace(1,0,256)))
cmaps['W2B'] = ListedColormap(W2B[0], name='W2B')

W2B_8 = np.dstack((np.linspace(1,0,8), np.linspace(1,0,8), np.linspace(1,0,8)))
cmaps['W2B_8'] = ListedColormap(W2B_8[0], name='W2B_8')

G2B = np.dstack((np.linspace(0,1,256), np.linspace(100/255,1,256), np.linspace(0,0,256)))
cmaps['G2B'] = ListedColormap(G2B[0], name='G2B')

R = np.concatenate((np.linspace(0,1,5), np.ones(10)))
G = np.concatenate((np.linspace(0,1,5), np.linspace(1,0,10)))
B = np.concatenate((np.ones(5), np.linspace(1,0,10)))
BWR = np.dstack((R, G, B))
cmaps['BWR'] = ListedColormap(BWR[0], name='BWR')
