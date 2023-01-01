import numpy as np
import equations as equ
from scipy.optimize import fsolve
import math
from  librerieTesi.diffuseRaman import core
import matplotlib.pyplot as plt
from  diffuse_optics import domodule
from  diffuse_optics import opticalGeomParam as o
import pandas as pd
from scipy.optimize import leastsq
v = 29.97925/1.4
pi = math.pi
exp =np.exp
log = np.log    
sqrt = np.sqrt
A =2.94852  
def derivative(x,y):
    dx = x[1]-x[0]
    return np.gradient(y, dx)

def de(rho,  z, t, z_e, mu_s):
    return (equ.diff_eq(t, rho, z, v, mu_s), equ.J_diff_eq(t, rho, z, v, mu_s), equ.d_de_dt(t, rho, z, v, mu_s))
def te(rho, z, t, z_e, mu_s):
    x = (t*v*mu_s)**2-3*(rho**2+z_e**2)*mu_s**2
    if np.isscalar(x) or len(x) == 1:
        if x<= 0:
            return (0,0,0)
        return (equ.te_c(t, rho, z, v, mu_s), equ.te_c_J(t, rho, z, v, mu_s), equ.d_te_dt_c(t, rho, z, v, mu_s))
    set_zero = (x<=0)
    res_phi = equ.te_c(t, rho, z, v, mu_s)
    res_J = equ.te_c_J(t, rho, z, v, mu_s)
    res_d_phi_dt = equ.d_te_dt_c(t, rho, z, v, mu_s)
    res_phi[set_zero] = 0  
    res_J[set_zero] = 0
    res_d_phi_dt[set_zero] = 0
    return (res_phi, res_J, res_d_phi_dt)
def te_M( rho, z, t, z_e, mu_s):
    x = (t*v*3*mu_s)**2-(rho**2+z_e**2)*(3*mu_s)**2
    if np.isscalar(x) or len(x) == 1:
        if x<= 0:
            return (0,0,0)
        return (equ.te_M(t, rho, z, v, mu_s), equ.te_M_J(t, rho, z, v, mu_s), equ.d_te_dt_M(t, rho, z, v, mu_s))
    set_zero = (x<=0)
    res_phi = equ.te_M(t, rho, z, v, mu_s)
    res_J = equ.te_M_J(t, rho, z, v, mu_s)
    res_d_phi_dt = equ.d_te_dt_M(t, rho, z, v, mu_s)
    res_phi[set_zero] = 0  
    res_J[set_zero] = 0
    res_d_phi_dt[set_zero] = 0
    return (res_phi, res_J, res_d_phi_dt)
def z_const_BC( phi, J,d_phi_dt, t, fit_fun, mu_s):
    return 2*A/(3*mu_s)
def pcbc_BC( phi, J,d_phi_dt, t, fit_fun, mu_s):
    res = []
    D = 1/(3*mu_s)
    z_e0 = 2*A*D

    for i in range(len(t)):
        find_ze = lambda z_e: fit_fun(phi( z_e, t[i])) - fit_fun(2*A*D*J(z_e, t[i]))
        res.append(fsolve(find_ze,z_e0*2)[0])
    z_e = np.array(res)
    return z_e
def pcbc_te_BC( phi, J,d_phi_dt, t, fit_fun, mu_s):
    res = []
    D = 1/(3*mu_s)
    z_e0 = 2*A*D
    for i in range(len(t)):
        find_ze = lambda z_e: fit_fun(phi( z_e, t[i])+d_phi_dt(z_e, t[i])/(v*mu_s)) - fit_fun(2*A*D*J(z_e, t[i]))
        res.append(fsolve(find_ze,z_e0*2)[0])
    z_e = np.array(res)
    return z_e
def PCBC(  phi, J, z_e, t, mu_s):
    return phi
def fick(  phi, J, z_e, t, mu_s):
    return J
def exp_PCBC( phi, J, z_e, t, mu_s):
    n_before = int(t[0]/(t[1]-t[0]))
    t_step = np.linspace(0,t[0], n_before)[:(n_before)]
    t_new = np.insert(t, 0, t_step)
    conv_mat = core.conv_matrix(np.exp(-mu_s*v*t_new), len(t_new))
    phi_new = np.insert(phi,0,[0]*n_before)
    res = np.matmul(conv_mat, phi_new)
    return res[n_before:]
    #return np.matmul(conv_mat, phi)
def exp_J(  phi, J, z_e, t, mu_s):
    n_before = int(t[0]/(t[1]-t[0]))
    t_step = np.linspace(0,t[0], n_before)[:(n_before)]
    t_new = np.insert(t, 0, t_step)
    conv_mat = core.conv_matrix(np.exp(-mu_s*v*t_new), len(t_new))
    J_new = np.insert(J,0,[0]*n_before)
    res = np.matmul(conv_mat, J_new)
    return res[n_before:]
def J_D(  phi, J, z_e, t, mu_s):
    return J- derivative(t,J)/(mu_s*v)
def PCBC_D(  phi, J, z_e, t, mu_s):
    return phi- derivative(t,phi)/(mu_s*v)
def J_D_plus(  phi, J, z_e, t, mu_s):
    return J+ derivative(t,J)/(mu_s*v)
def PCBC_D_plus(  phi, J, z_e, t, mu_s):
    return phi+ derivative(t,phi)/(mu_s*v)
        

def forward(rho, time, eq, bc, f, mu_s,fit_fun):
    phi = lambda  z_e, time : eq(rho = rho, t = time, z = 0, z_e = z_e, mu_s = mu_s)[0] - eq(rho = rho, t = time, z = z_e, z_e = z_e, mu_s = mu_s)[0]#TODO forse la H deve essere qui
    J  = lambda  z_e, time : eq(rho = rho, t = time, z_e = z_e, z = z_e, mu_s = mu_s)[1] #usare flux come PCBC!! perchè ci sono combinazioni che si rieptono
    d_phi_dt = lambda  z_e, time : eq(rho = rho, t = time, z = 0, z_e = z_e, mu_s = mu_s)[2] - eq(rho = rho, t = time, z = z_e, z_e = z_e, mu_s = mu_s)[2]
    z_e = bc(phi, J, d_phi_dt, time,fit_fun, mu_s = mu_s)
    tof = f(phi(z_e, time), J(z_e, time), z_e, time, mu_s = mu_s)
    t0 = rho/v
    #tof[time<t0] = 0 #TODO sperimentale
    tof[tof<0] = 0
    return tof#/np.sum(tof)

