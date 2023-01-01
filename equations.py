import numpy as np
from scipy.special import iv
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
pi = math.pi
exp =np.exp
sqrt = np.sqrt
#TODO togliere  z_e = 2*z_e per te perchè instabile & log nel main
def te_c(t, rho, z_e, v, mu_s):
    #C = (mu_s/3)**(3/2)/(8*math.pi*v)
    #z_e = 2*z_e
    C = (v/(8*pi))*(3*mu_s**2)**(3/2)
    x = (t*v*mu_s)**2-3*(rho**2+z_e**2)*mu_s**2
    
    if np.isscalar(x):
        if x<=0:
            return 0     
        x = sqrt(x)
        if (x/2)> 100:
            return te_1(t, rho, z_e, v, mu_s)
        res = C*iv(1,x/2) * np.exp(-t*v*mu_s/2)/x
        return res
    set_zero = x<=0
    x = sqrt(x)
    res = C*iv(1,x/2) * np.exp(-t*v*mu_s/2)/x
    res[set_zero] = 0
    use_exp = (x/2)> 100
    res_exp = te_1(t, rho, z_e, v, mu_s)
    res[use_exp] = res_exp[use_exp]
    return res
def te_1(t,rho, z_e, v, mu_s):
    #z_e = 2*z_e
    C = (v/(8*pi))*(3*mu_s**2)**(3/2)
    x = (t*v*mu_s)**2-3*(rho**2+z_e**2)*mu_s**2 
    x = sqrt(x)
    mbf = lambda x , x2 : exp(x-x2)/sqrt(2*pi*x)
    return  C*mbf(x/2, t*v*mu_s/2)/x
def te_c_J(t, rho, z_e, v, mu_s):
    #z_e = 2*z_e
    x = (t*v*mu_s)**2-3*(rho**2+z_e**2)*mu_s**2
    C = (v/(8*pi))*(3*mu_s**2)**(3/2)
    if np.isscalar(x):
        if x<=0:
            return 0  
        x = sqrt(x)
        if (x/2)> 100:
            return - te_c_J_1(t, rho, z_e, v, mu_s)
        C = (v/(8*pi))*(3*mu_s**2)**(3/2)
        res = C*exp(-v*mu_s*t/2)*3*mu_s**2*z_e*(iv(1,x/2)/x**3 -(iv(0,x/2) + iv(2,x/2))/(4*x**2))
        return -res
    set_zero = (x<=0)
    x = sqrt(x)
    res = C*exp(-v*mu_s*t/2)*3*mu_s**2*z_e*(iv(1,x/2)/x**3 -(iv(0,x/2) + iv(2,x/2))/(4*x**2))
    res[set_zero] = 0   
    use_exp = (x/2)> 100
    res_exp = te_c_J_1(t, rho, z_e, v, mu_s)
    res[use_exp] = res_exp[use_exp]
    return -res
def te_c_J_1(t, rho, z_e, v, mu_s):
    #z_e = 2*z_e#TODO ERRORI, per ora non correggo perchè sembra stabile
    x = (t*v*mu_s)**2-3*(rho**2+z_e**2)*mu_s**2
    x = sqrt(x)
    C = (v/(8*pi))*(3*mu_s**2)**(3/2)
    mbf = lambda x , x2 : exp(x-x2)/sqrt(2*pi*x)
    res = C*3*mu_s**2*z_e*(1/x**3 -1/(2*x**2))*mbf(x/2,v*mu_s*t/2)
    return res
def d_te_dt_c(t, rho, z_e, v, mu_s):
    C = (v/(8*pi))*(3*mu_s**2)**(3/2)
    a = (mu_s*v)
    b = 3*(rho**2+z_e**2)*mu_s**2
    x = (a*t)**2- b
    if np.isscalar(x):
        if x<=0:
            return 0     
        if (x/2)> 100:
            return d_te_dt_c_1(t, rho, z_e, v, mu_s)
        P1 = - a**2 *t * iv(1, sqrt(x)/2)/x**(3/2)
        P2 = a**2*t *(iv(0,sqrt(x)/2)+iv(2,sqrt(x)/2))/(4*x)
        P3 = - a*iv(1,sqrt(x)/2)/(2*sqrt(x))
        res = C*exp(-a*t/2)*(P1+P2+P3)
        return res
    P1 = - a**2 *t * iv(1, sqrt(x)/2)/x**(3/2)
    P2 = a**2*t *(iv(0,sqrt(x)/2)+iv(2,sqrt(x)/2))/(4*x)
    P3 = - a*iv(1,sqrt(x)/2)/(2*sqrt(x))
    res = C*exp(-a*t/2)*(P1+P2+P3)
    set_zero = x<=0
    res[set_zero] = 0
    use_exp = (x/2)> 100
    res_exp = d_te_dt_c_1(t, rho, z_e, v, mu_s)
    res[use_exp] = res_exp[use_exp]
    return res
def d_te_dt_c_1(t, rho, z_e, v, mu_s):
    C = (v/(8*pi))*(3*mu_s**2)**(3/2)
    a = (mu_s*v)
    b = 3*(rho**2+z_e**2)*mu_s**2
    x = (a*t)**2- b
    
    
    res = exp(sqrt(x)/2-a*t/2)*((a**2*t/(2*sqrt(x)) - a/2)/x**(3/4) - 3*a**2*t/(2*x**(7/4)))/sqrt(pi)
    return C*res
    
    
    
    
def te_M(t, rho, z_e, v, mu_s):
    #C = (mu_s/3)**(3/2)/(8*math.pi*v)
    #z_e = 2*z_e
    C = v*(3*mu_s)**3/(8*pi)
    
    x = (t*v*3*mu_s)**2-(rho**2+z_e**2)*(3*mu_s)**2
    #plt.plot(sqrt(x)/2-t*v*(3*mu_s)/2)
    if np.isscalar(x):
        if x<=0:
            return 0     
        x = sqrt(x)
        if (x/2)> 100    :
            return te_M_1(t, rho, z_e, v, mu_s)
        res = C*iv(1,x/2) * np.exp(-t*v*3*mu_s/2)/x
        return res
    set_zero = x<=0
    x = sqrt(x)
    res = C*iv(1,x/2) * np.exp(-t*v*3*mu_s/2)/x
    res[set_zero] = 0
    use_exp = (x/2)> 100
    res_exp = te_M_1(t, rho, z_e, v, mu_s)
    res[use_exp] = res_exp[use_exp]
    return res
def te_M_1(t,rho, z_e, v, mu_s):
    #z_e = 2*z_e
    C = v*(3*mu_s)**3/(8*pi**(3/2))
    x = (t*v*3*mu_s)**2-(rho**2+z_e**2)*(3*mu_s)**2 
    #mbf = lambda x , x2 : exp(x-x2)/sqrt(2*pi*x)
    return  C*exp(sqrt(x)/2-(3*mu_s*v*t)/2)/x**(3/4)#C* mbf(x/2, t*v*3*mu_s/2)/x
def te_M_J(t, rho, z_e, v, mu_s):
    #z_e = 2*z_e
    C = v*(3*mu_s)**3/(8*pi)
    a = (3*mu_s*v)**2 - (3*rho*mu_s)**2
    b = (3*mu_s)**2
    x = a-b*z_e**2
    
    if np.isscalar(x):
        if x<=0:
            return 0  
        #x = sqrt(x)
        if (x/2)> 100:
            return te_M_J_1(t, rho, z_e, v, mu_s)
        res = C*exp(-v*(3*mu_s)*t/2)*b*z_e*(iv(1,sqrt(x)/2)/x**(3/2) - (iv(0,sqrt(x)/2) +iv(2,sqrt(x)/2))/(4*x))
        return res
    set_zero = (x<=0)
    #x = sqrt(x)
    res = C*exp(-v*(3*mu_s)*t/2)*b*z_e*(iv(1,sqrt(x)/2)/x**(3/2) - (iv(0,sqrt(x)/2) +iv(2,sqrt(x)/2))/(4*x))
    res[set_zero] = 0   
    use_exp = (x/2)> 100
    res_exp = te_M_J_1(t, rho, z_e, v, mu_s)
    res[use_exp] = res_exp[use_exp]
    return res
def te_M_J_1(t, rho, z_e, v, mu_s):
    #z_e = 2*z_e#TODO ERRORI, per ora non correggo perchè sembra stabile
    C = v*(3*mu_s)**3/(8*pi**(3/2))
    a = (3*mu_s*v*t)**2 - (3*rho*mu_s)**2
    b = (3*mu_s)**2
    x = a-b*z_e**2
    
    mbf = lambda x , x2 : exp(x-x2)/sqrt(2*pi*x)
    res = C*(3*b*z_e/(2*x**(7/4)) -b*z_e/(2*x**(5/4)))*mbf(sqrt(x)/2,v*(3*mu_s)*t/2)
    res = C*(-b*z_e/(2*x**(5/4)))*exp(sqrt(x)/2-(3*mu_s*v*t)/2)
    return -res
def d_te_dt_M(t, rho, z_e, v, mu_s):
    C = v*(3*mu_s)**3/(8*pi**(3/2))
    a = (3*mu_s*v)
    b = (rho**2+z_e**2)*(3*mu_s)**2
    x = (a*t)**2- b
    if np.isscalar(x):
        if x<=0:
            return 0     
        if (x/2)> 100:
            return d_te_dt_1_M(t, rho, z_e, v, mu_s)
        P1 = - a**2 *t * iv(1, sqrt(x)/2)/x**(3/2)
        P2 = a**2*t *(iv(0,sqrt(x)/2)+iv(2,sqrt(x)/2))/(4*x)
        P3 = - a*iv(1,sqrt(x)/2)/(2*sqrt(x))
        res = exp(-a*t/2)*(P1+P2+P3)
        return C*res
    P1 = - a**2 *t * iv(1, sqrt(x)/2)/x**(3/2)
    P2 = a**2*t *(iv(0,sqrt(x)/2)+iv(2,sqrt(x)/2))/(4*x)
    P3 = - a*iv(1,sqrt(x)/2)/(2*sqrt(x))
    res = C*exp(-a*t/2)*(P1+P2+P3)
    set_zero = x<=0
    res[set_zero] = 0
    use_exp = (x/2)> 100
    res_exp = d_te_dt_1_M(t, rho, z_e, v, mu_s)
    res[use_exp] = res_exp[use_exp]
    return res
def d_te_dt_1_M(t, rho, z_e, v, mu_s):
    
    a = (3*mu_s*v)
    b = (rho**2+z_e**2)*(3*mu_s)**2
    x = (a*t)**2- b
    C = v*(3*mu_s)**3/(8*pi**(3/2))
    
    res = exp(sqrt(x)/2-a*t/2)*((a**2*t/(2*sqrt(x)) - a/2)/x**(3/4) - 3*a**2*t/(2*x**(7/4)))/sqrt(pi)
    return C*res
def diff_eq(t, rho, z_e, v, mu_s):
    D = 1/(3*mu_s)
    #plt.plot((rho**2+(2*z_e)**2)/(4*D*v*t))
    
    res = v*exp(-(rho**2+(2*z_e)**2)/(4*D*v*t))/(4*pi*D*v*t)**(3/2)
    #res = v*exp(-(4*z_e**2)/(4*D*v*t))
    return res

def J_diff_eq(t, rho, z_e, v, mu_s):
    D = 1/(3*mu_s)#controllare z_e
    res = ((2*z_e)/(2*D*t))*exp(-(rho**2+(2*z_e)**2)/(4*D*v*t))/(4*pi*D*v*t)**(3/2)
    #res = (z_e/(D*t))*exp(-(4*z_e**2)/(4*D*v*t))
    #res = ((2*z_e)/(2*D*t))/(4*pi*D*v*t)**(3/2)
    return res
def d_de_dt(t, rho, z_e, v, mu_s):
    D = 1/(3*mu_s)
    C = v/(4*pi*D*v)**(3/2)
    a = (mu_s*v)
    
    b = 3*(rho**2+(2*z_e)**2)*mu_s**2
    res =  exp(-b/(4*a*t))*(b-6*a*t)/(4*a*t**(7/2))
    #res = v*exp(-(4*z_e**2)/(4*D*v*t))
    return C*res
if __name__ == "__main__":
    mu_s = 10
    v = 29.97925/1.4

    z_e = 0.1966
    A =2.94852
    rho = 1
    D = 1/(3*mu_s)
    def phi(t, rho, z):
    
        x = (t*v*mu_s)**2-3*((rho*mu_s)**2+(z*mu_s)**2)
        x = sqrt(x)
        res = te_M(t,rho,z, v, mu_s)
        return res
    
    t = np.arange(5,4000)*1e-3
    
    res =  te_M(t, rho, 2*z_e, v, mu_s)
    plt.plot(res)
    res = diff_eq(t, rho, z_e, v, mu_s)
    plt.plot(res)
    plt.show()

