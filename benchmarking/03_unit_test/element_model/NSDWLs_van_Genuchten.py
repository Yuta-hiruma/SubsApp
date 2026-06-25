# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 13:03:44 2022

@author: tomdegast
"""
import numpy as np
import scipy as sci

def van_genuchten (WRC, theta):
    #iterertavely slove psi
    Se = np.zeros(np.shape(theta))
    psi_out = np.zeros(np.shape(theta))
    #shrinking
    Se = (theta - WRC.vG_WCR) / (WRC.vG_WCS - WRC.vG_WCR)
    #print (Se)
    #psi_out = psi_in
    #print (psi_out)
    for i in range(len(Se)):     
        psi_high = WRC.vG_a[i]
        psi_low = 0.0
        psi_avg = 0.5*psi_high+0.5*psi_low
        tolerance = np.abs(psi_high-psi_low)
        if Se[i] >= 1.0: 
            Se[i] = 1.0
            psi_out[i] = 0.0
            tolerance = 0.0
        if Se[i] <= 0.0 :
            Se[i] = 0.0
            psi_out[i] = WRC.vG_a[i]
            tolerance = 0.0
        while tolerance > 0.001:
            Se_high = (1-(np.log(1+psi_high/WRC.vG_a[i]))/np.log(2))*((1/(1+(WRC.vG_alpha[i]*psi_high)**WRC.vG_n[i]))**WRC.vG_m[i])
            Se_avg = (1-(np.log(1+psi_avg/WRC.vG_a[i]))/np.log(2))*((1/(1+(WRC.vG_alpha[i]*psi_avg)**WRC.vG_n[i]))**WRC.vG_m[i])
            Se_low = (1-(np.log(1+psi_low/WRC.vG_a[i]))/np.log(2))*((1/(1+(WRC.vG_alpha[i]*psi_low)**WRC.vG_n[i]))**WRC.vG_m[i])
            if Se_high == Se[i] : 
                psi_out[i] = psi_high 
                tolerance = 0.0
            if Se_avg == Se[i] : 
                psi_out[i] = psi_avg 
                tolerance = 0.0
            if Se_low == Se[i] : 
                psi_out[i] = psi_low 
                tolerance = 0.0
            psi_out[i] = psi_avg
            tolerance = np.abs(psi_high-psi_low)
            #print (i, tolerance,"psi", psi_low, psi_avg, psi_high,"Se_avg", Se_avg,"Se[i]",Se[i])
            if Se_avg < Se[i] : 
                psi_high = psi_avg
                psi_avg = 0.5*psi_high+0.5*psi_low
            if Se_avg > Se[i] : 
                psi_low = psi_avg
                psi_avg = 0.5*psi_high+0.5*psi_low
            
        #print(i,Se[i],psi_out[i])
           
    #print(psi_out)
    return psi_out    

def derrivative_van_genuchten_2 (WRC, theta, psi_m):
    Se = np.zeros(np.shape(theta))
    dtheta_dpsi = np.zeros(np.shape(theta))
    psi_m = -psi_m
    #shrinking
    
    for i in range(len(dtheta_dpsi)):
        if psi_m[i] <= 0 :
            psi_b = psi_m[i]
            psi_c = psi_m[i] +0.005
            Se_b = (1-(np.log(1+psi_b/WRC.vG_a[i]))/np.log(2))*((1/(1+(WRC.vG_alpha[i]*psi_b)**WRC.vG_n[i]))**WRC.vG_m[i])
            Se_c = (1-(np.log(1+psi_c/WRC.vG_a[i]))/np.log(2))*((1/(1+(WRC.vG_alpha[i]*psi_c)**WRC.vG_n[i]))**WRC.vG_m[i])
            theta_b = (WRC.vG_WCS[i]-WRC.vG_WCR[i])*Se_b+WRC.vG_WCR[i]
            theta_c = (WRC.vG_WCS[i]-WRC.vG_WCR[i])*Se_c+WRC.vG_WCR[i]
            
            dtheta_dpsi[i] = (theta_c-theta_b)/(2*0.005)
        else :
            psi_a = psi_m[i] -psi_m[i]*0.5
            psi_c = psi_m[i] +psi_m[i]*0.5
            Se_a = (1-(np.log(1+psi_a/WRC.vG_a[i]))/np.log(2))*((1/(1+(WRC.vG_alpha[i]*psi_a)**WRC.vG_n[i]))**WRC.vG_m[i])
            Se_c = (1-(np.log(1+psi_c/WRC.vG_a[i]))/np.log(2))*((1/(1+(WRC.vG_alpha[i]*psi_c)**WRC.vG_n[i]))**WRC.vG_m[i])
            theta_a = (WRC.vG_WCS[i]-WRC.vG_WCR[i])*Se_a+WRC.vG_WCR[i]
            theta_c = (WRC.vG_WCS[i]-WRC.vG_WCR[i])*Se_c+WRC.vG_WCR[i]
        
            dtheta_dpsi[i] = (theta_c-theta_a)/(2*psi_m[i]*0.5)
   
   
    
    #print("psi_m", psi_m)
    
    return dtheta_dpsi


def derrivative_van_genuchten (WRC, theta,e, psi_m):
    Se = np.zeros(np.shape(theta))
    dSe_dpsi = np.zeros(np.shape(theta))
    dtheta_dpsi = np.zeros(np.shape(theta))
    psi_m = -psi_m
    #shrinking
    a = WRC.vG_a
    alpha = WRC.vG_alpha
    n = WRC.vG_n
    m = WRC.vG_m
        
    dSe_dpsi =(alpha*(-m)*n*(1-(np.log((psi_m/a)+1)/np.log(2))))*((alpha*psi_m)**(n-1))*(((alpha*psi_m)**n+1)**(-m-1))-((((alpha*psi_m)**n+1)**(-m))/(a*np.log(2)*((psi_m/a)+1)))
            
    #dtheta_dpsi = dSe_dpsi*(WRC.vG_WCS-WRC.vG_WCR)
    dtheta_dpsi = dSe_dpsi*e
    
    for i in range(len(psi_m)):
        if psi_m[i] <= 0:
            dtheta_dpsi[i] = 0.0
    #print("psi_m", psi_m)
    
    return dtheta_dpsi

def derrivative_van_genuchten_S_e (WRC, theta,e, psi_m):
    Se = np.zeros(np.shape(theta))
    dSe_dpsi = np.zeros(np.shape(theta))
    dtheta_dpsi = np.zeros(np.shape(theta))
    psi_m = -psi_m
    #shrinking
    a = WRC.vG_a
    alpha = WRC.vG_alpha
    n = WRC.vG_n
    m = WRC.vG_m
        
    dSe_dpsi =(alpha*(-m)*n*(1-(np.log((psi_m/a)+1)/np.log(2))))*((alpha*psi_m)**(n-1))*(((alpha*psi_m)**n+1)**(-m-1))-((((alpha*psi_m)**n+1)**(-m))/(a*np.log(2)*((psi_m/a)+1)))
            
    #dtheta_dpsi = dSe_dpsi*(WRC.vG_WCS-WRC.vG_WCR)
    #dtheta_dpsi = dSe_dpsi*e
    
    for i in range(len(psi_m)):
        if psi_m[i] <= 0:
            dSe_dpsi[i] = 0.0
    #print("psi_m", psi_m)
    
    return dSe_dpsi

def derrivative_van_genuchten_S_e_2(WRC, theta,e, psi_m):
    Se = np.zeros(np.shape(theta))
    dSe_dpsi = np.zeros(np.shape(theta))
    dtheta_dpsi = np.zeros(np.shape(theta))
    psi_m = -psi_m
    #shrinking
    a = WRC.vG_a
    alpha = WRC.vG_alpha
    n = WRC.vG_n
    m = WRC.vG_m
        
    dSe_dpsi =(alpha*(-m)*n*(1-(np.log((psi_m/a)+1)/np.log(2))))*((alpha*psi_m)**(n-1))*(((alpha*psi_m)**n+1)**(-m-1))-((((alpha*psi_m)**n+1)**(-m))/(a*np.log(2)*((psi_m/a)+1)))
            
    #dtheta_dpsi = dSe_dpsi*(WRC.vG_WCS-WRC.vG_WCR)
    #dtheta_dpsi = dSe_dpsi*e
    
    for i in range(len(psi_m)):
        if psi_m[i] <= 0 : dSe_dpsi[i] = 0.0

    #print("psi_m", psi_m)
    
    return dSe_dpsi


def derrivative_van_genuchten_S_e_loc(WRC, theta,e, psi_m):
    Se = 0
    dSe_dpsi = 0
    dtheta_dpsi = 0
    psi_m = -psi_m
    #shrinking
    a = WRC.vG_a
    alpha = WRC.vG_alpha
    n = WRC.vG_n
    m = WRC.vG_m
        
    dSe_dpsi =(alpha*(-m)*n*(1-(np.log((psi_m/a)+1)/np.log(2))))*((alpha*psi_m)**(n-1))*(((alpha*psi_m)**n+1)**(-m-1))-((((alpha*psi_m)**n+1)**(-m))/(a*np.log(2)*((psi_m/a)+1)))
            
    #dtheta_dpsi = dSe_dpsi*(WRC.vG_WCS-WRC.vG_WCR)
    #dtheta_dpsi = dSe_dpsi*e
    
    if psi_m <= 0 : Se_dpsi = 0.0

    #print("psi_m", psi_m)
    
    return dSe_dpsi


"""
    psi_n = ((1/((1/Se)**WRC.vG_m))**WRC.vG_n)/WRC.vG_alpha
    
    for i in range(len(Se)):
        rel_tol =1
        Cs_loop = 0
        Se_loop = 0
        j=1
        psi_old = psi_n[i]
        while rel_tol > 0.01:
            if psi_old >WRC.vG_a[i] : psi_old = WRC.Vg_a[i]
            Cs_loop = (1+(np.log(1+psi_old/WRC.vG_a[i]))/np.log(2))
            Se_loop = Cs_loop*(1/ (1+(WRC.vG_alpha[i]*psi_old)**WRC.vG_n[i]))**WRC.vG_m[i]
            Cs_loop_low = (1+(np.log(1+0.99*psi_old/WRC.vG_a[i]))/np.log(2))
            Se_loop_low = Cs_loop_low*(1/ (1+(WRC.vG_alpha[i]*psi_old*.99)**WRC.vG_n[i]))**WRC.vG_m[i]
            ds_dSe = (psi_old*.99-psi_old)/(Se_loop_low-Se_loop)
            #if ds_dSe < 0 : ds_dSe=0
            dSe = Se[i] - Se_loop
            ds = ds_dSe *dSe
            psi_old = psi_old + ds*0.5
            print (psi_old)
            rel_tol = (psi_n[i]-psi_old)/psi_n[i]
            psi_n[i] = psi_old
"""  
