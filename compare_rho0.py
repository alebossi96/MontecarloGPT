import numpy as np
import equations as equ
from scipy.optimize import fsolve
import math
from  librerieTesi.diffuseRaman import core
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import core_functions as cf
from scipy.special import iv
import montecarlomodule as mc
v = 29.97925#/1.4
pi = math.pi
exp =np.exp
log = np.log    
sqrt = np.sqrt
def diff_eq(t, r, mu_s):
    D = 1/(3*mu_s)
    res = v/(4*pi*D*v*t)**(3/2)*exp(-r**2/(4*D*v*t))
    if t[0] == 0:
        res[0] = 0
    return res
def diff_eq_corrEXP(t, r, mu_s):
    de = diff_eq(t, r, mu_s)
    conv_mat = core.conv_matrix(np.exp(-mu_s*v*t), len(t))
    return np.matmul(conv_mat, de)
def derivative_t(v, D, t, r):
    return (np.exp(-r**2/(4*D*t*v)) * (r**2 - 6*D*t*v))/(4*D*t**(7/2)*v)
def diff_eq_corrDER(t, r, mu_s):
    de = diff_eq(t, r, mu_s)
    D = 1/(3*mu_s)
    d_de_dt = derivative_t(v, D, t, r)
    res =  de - 1/(mu_s*v)*d_de_dt
    res[res<0] = 0
    return res
class obj_det:
    def __init__(self, x,y,z, r):
        self.x = x
        self.y = y
        self.z = z
        self.r = r
d = obj_det(1,0,0,0.05)
g = 0.5
for mu_s in [2]:#[1,4,8,16,27,81]:
    res = mc.mc(g,mu_s, d)
    time = np.arange(0,len(res))/1e3
    plt.plot(time, res/np.sum(res[100:]), label = "$\mu_{s,mc} = $" +str(mu_s)+ "g = "+ str(g))
    res = diff_eq_corrEXP(time, 1, mu_s*(1-g))
    plt.plot(time, res/np.sum(res[100:]), label = "$\mu_{s,exp} =$" +str(mu_s)+" g "+ str(g))
    res = diff_eq(time, 1, mu_s*(1-g))
    plt.plot(time, res/np.sum(res[100:]), label = "$\mu_{s,de} = $" +str(mu_s)+" g "+ str(g))
    #res = diff_eq_corrDER(time, 0, mu_s*(1-g))
    #plt.plot(time, res/np.sum(res[100:]), label = "$\mu_{s,der} = $" +str(mu_s)+" g "+ str(g))
    #
    #
    #plt.savefig("imgs/"+str(mu_s)+".png")
plt.yscale('log')
plt.legend()
plt.show()
#    plt.close()
#    plt.show()
    
