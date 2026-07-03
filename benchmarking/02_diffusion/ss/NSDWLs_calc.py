
""" 
@author: tomdegast

210723 : This is the first attempt to simulate the paper written by Phil et al in Python
210902 : Received the Matlabscript by Phil et al - rewriting to Python

"""


#  Import numpy and the ODE solver from Scipy
import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import odeint

from NSDWLs_equations import *

#  Import plotting functionality
import matplotlib.pyplot as plt


class out:
    def __init__ (self):  
        self.t = []
        self.y = []


def RK4_mine (equations,simulation_time,input_in,method,t_eval,max_step, rtol):
    
    output = out()
    output.t = t_eval
    output.y = np.zeros([len(input_in),len(t_eval)])
    
    h = t_eval[1]-t_eval[0]
    print(h)
    while h > max_step : h=h/2
    print(h)
    
    simssteps = int(t_eval[-1]/h)
    
    t_temp = np.linspace (t_eval[0],t_eval[-1],simssteps+1) 
    y_temp = np.zeros([len(input_in),len(t_temp)])
    
    output.y[:,0] = input_in
    y_temp[:,0] = input_in
    j=0
    
    
    for i in range(len(t_temp)-1):
        #print(y_temp[:,i])
        k1 = h * np.array(equations(t_temp[i], y_temp[:,i].copy()))
        k2 = h * np.array(equations(t_temp[i] + h/2, y_temp[:,i].copy() + k1/2))
        k3 = h * np.array(equations(t_temp[i] + h/2, y_temp[:,i].copy() + k2/2))
        k4 = h * np.array(equations(t_temp[i] + h, y_temp[:,i].copy() + k3))
        
        # Update the next value of y and t
        
        y_temp[:,i+1] = y_temp[:,i] + (k1 + 2*k2 + 2*k3 + k4) / 6


        #print(k1[0,40], k1[0,80])
        #print(k2[0,40], k2[0,80])
        #print(k3[0,40], k3[0,80])
        #print(k4[0,40], k4[0,80])
        #print(y_temp[40,i], y_temp[80,i])
        #print(y_temp[40,i+1], y_temp[80,i+1])

        if t_temp[i+1] == t_eval[j+1] : 
            output.y[:,j+1] = y_temp[:,i+1] 
            j = j+1    
         
    
    return output

def RK8_mine(equations,simulation_time,input_in,method,t_eval,max_step, rtol):
    output = out()
    output.t = t_eval
    output.y = np.zeros([len(input_in),len(t_eval)])
    
    h = t_eval[1]-t_eval[0]
    print(h)
    while h > max_step : h=h/2
    print(h)
    
    simssteps = int(t_eval[-1]/h)
    
    t_temp = np.linspace (t_eval[0],t_eval[-1],simssteps+1) 
    y_temp = np.zeros([len(input_in),len(t_temp)])
    
    output.y[:,0] = input_in
    y_temp[:,0] = input_in
    j=0  
    
    # Butcher tableau coefficients for RK8 method
    # Coefficients come from the classical RK8 method
    c = np.array([0, 1/2, 1/2, 1, 1/6, 1/3, 1/3, 1/2, 1/6, 1/6, 1/3, 1/3, 1/2])
    
    a = np.array([
        [],  # c[0] = 0
        [1/2],  # c[1] = 1/2
        [0, 1/2],  # c[2] = 1/2
        [0, 0, 1],  # c[3] = 1
        [1/6, 1/6, 0, 0],  # c[4] = 1/6
        [1/3, -1/3, 1/3, 0],  # c[5] = 1/3
        [1/3, 0, 1/6, -1/6, 1/6],  # c[6] = 1/3
        [1/2, 1/6, 0, 1/6, 0, 1/6],  # c[7] = 1/2
        [1/6, 1/3, 0, 1/6, 0, 1/3],  # c[8] = 1/6
        [1/6, 0, 0, 1/3, 0, 1/3],  # c[9] = 1/6
        [1/3, 0, 1/6, -1/6, 0, 1/3],  # c[10] = 1/3
        [1/3, 1/3, 0, -1/6, 0, 1/3],  # c[11] = 1/3
        [1/2, 0, 1/6, 0, -1/6, 1/3],  # c[12] = 1/2
    ])
    
    b = np.array([1/90, 0, 4/90, 4/90, 4/90, 4/90, 4/90, 4/90, 1/90, 0])

    t_values = np.arange(t0, t_end + h, h)  # Time grid
    y_values = [y0]  # Initial value of y

    y_current = y0
    for t in t_values[:-1]:
        k = []
        k.append(h * f(t, y_current))  # k1
        for i in range(1, len(c)):
            t_stage = t + c[i] * h
            y_stage = y_current + sum(a[i][j] * k[j] for j in range(i))
            k.append(h * f(t_stage, y_stage))

        # Update y using the weighted average of all the stages
        y_next = y_current + sum(b[i] * k[i] for i in range(len(b)))
        y_values.append(y_next)
        
        # Update the current value of y for the next iteration
        y_current = y_next

    return t_values, np.array(y_values)


def equations_test4(sim_time, input_in):
    output = np.zeros(np.shape(input_in))
    
    output = np.cos(output+sim_time)
    
    return output


#simulation_time = [0, 100]
#input_in = [0,0]

#method = 'RK_mine'
#simssteps = 1000
#t_eval = np.linspace (simulation_time[0],simulation_time[-1],simssteps+1) 
#max_step = 0.1
#rtol = 1e-12

#sol = RK_mine (equations_test4,simulation_time,input_in,method,t_eval,max_step, rtol)

#plt.figure()
#plt.plot(sol.t, sol.y[1,:])
#plt.plot(sol.t, np.sin(sol.t),'ok' )




