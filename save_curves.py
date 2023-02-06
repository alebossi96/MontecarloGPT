import matplotlib.pyplot as plt
import numpy as np
import math
import h5py
hf = h5py.File('data.h5', 'w')
class obj_det:
    def __init__(self, x,y,z, r):
        self.x = x
        self.y = y
        self.z = z
        self.r = r
import montecarlomodule as mc
mu_s_list = [1,5, 25, 125]
g_list = [0,0.5,0.9]
rho_list = [0,1]

for mu_s in mu_s_list:
    hf.create_group(str(mu_s))
    folder = hf.get(str(mu_s))
    for g in g_list:
        folder.create_group(str(g))

for mu_s in mu_s_list:
    fold1 =  hf [str(mu_s)]
    for g in g_list:
        fold2 = fold1[str(g)]
        for rho in rho_list:
            print(mu_s, g, rho)
            
            if rho>0.2:
                r = 0.1
            else:
                r = 0.01
            d = obj_det(rho,0,0,r)
            res = mc.mc(g,mu_s, d)
            fold2.create_dataset(str(rho),data=res)
hf.close()
