# Compressible Aerodynamic Calculator by Gabe Colangelo

import tkinter as tk
from tkinter import simpledialog
from tkinter import messagebox
import numpy as np
import math

def Isentropic(M,gamma):
    T_ratio = str((1 + ((gamma-1)/2)*M**2))
    P_ratio = str(((1 + ((gamma-1)/2)*M**2)**(gamma/(gamma-1))))
    rho_ratio = str(((1 + ((gamma-1)/2)*M**2)**(1/(gamma-1))))
    Mach_ang = str((np.arcsin(1/M))*(180/np.pi))
    Area_ratio = str(math.sqrt((1/M**2)*(((2+(gamma-1)*M**2)/(gamma+1))**((gamma+1)/(gamma-1)))))
    PM_angle = str((math.sqrt((gamma+1)/(gamma-1))*(np.arctan(math.sqrt((gamma-1)*(M**2-1)/(gamma+1))))-(np.arctan(math.sqrt(M**2-1))))*(180/np.pi))
    
    info = messagebox.showinfo('Isentropic Flow Results', 'For a Mach number of ' + str(M) + ' and a specific heat ratio of ' + str(gamma) + ':\n' +'\n' + 'T0/T = ' + T_ratio +'\n'+ 'P0/P = ' + P_ratio + '\n' + 'Mach Angle (deg) = ' + Mach_ang + '\n' + 'A/A* = ' + Area_ratio + '\n' + 'PM angle (deg) = ' + PM_angle, parent=ROOT)

def Normal(M,gamma):
    M2 = str(math.sqrt((1+((gamma-1)/2)*M**2)/(gamma*M**2 - (gamma-1)/2)))
    P2_P1 = str(1+(2*gamma*(M**2 -1 )/(gamma+1)))
    T2_T1 = str((1 + ((2*gamma*(M**2 -1))/(gamma+1)))*((2+(gamma-1)*M**2)/((gamma+1)*M**2)))
    P02_Pinf = str(((((gamma+1)*M**2)/2)**(gamma/(gamma-1)))/((2*gamma*M**2)/(gamma+1) - ((gamma-1)/(gamma+1)))**(1/(gamma-1)))
    
    info = messagebox.showinfo('Normal Shock Results', 'For a Mach number of ' + str(M) + ' and a specific heat ratio of ' + str(gamma) + ':\n' + '\n' + 'M2 = ' + M2 +'\n'+ 'P2/P1 = ' + P2_P1 + '\n' + 'T2/T1 = ' + T2_T1 + '\n' + 'P02/P1 = ' + P02_Pinf , parent=ROOT)

def Oblique(M,gamma,theta):
    n=2 # Weak Shock Solution
    theta_r = theta*np.pi/180
    a = (1+((gamma-1)/2)*M**2)*np.tan(theta_r)
    b = 1 - M**2
    c = (1+((gamma+1)/2)*M**2)*np.tan(theta_r)
    lamda = math.sqrt(b**2 - (3*a*c))
    chi = (9*a*b*c - (2*b**3) -(27*a**2))/(2*lamda**3)
    B = np.arctan((-b + 2*lamda*np.cos((np.arccos(chi)/3 + (2*np.pi*n/3))))/(3*a))
    beta = B*180/np.pi
    Mn1 = M*np.sin(B)
    Mn2 = math.sqrt((1+((gamma-1)/2)*Mn1**2)/(gamma*Mn1**2 -((gamma-1)/2)))
    M2 = str(Mn2/(np.sin( B-theta_r )))
    P2_P1 = str(1+(2*gamma*(Mn1**2 -1 )/(gamma+1)))
    T2_T1 = str((1 + ((2*gamma*(Mn1**2 -1))/(gamma+1)))*((2+(gamma-1)*Mn1**2)/((gamma+1)*Mn1**2)))
    
    info = messagebox.showinfo('Oblique Shock Results', 'For a Mach number of ' + str(M) + ', a specific heat ratio of ' + str(gamma) + ', and a weak turn angle of ' + str(theta) + ' degrees :\n' + '\n' + 'Shock Angle = '+ str(beta) +' degrees' + '\n'+ 'Mn1 = '+ str(Mn1) + '\n' + 'Mn2 = '+ str(Mn2) + '\n' + 'M2 = ' + M2 +'\n'+ 'P2/P1 = ' + P2_P1 + '\n' + 'T2/T1 = ' + T2_T1 , parent=ROOT)

def Rayleigh(M,gamma):
    P_PP = (1+gamma)/(1+(gamma*M**2))
    T_TT = M**2*((1+gamma)/(1+(gamma*M**2)))**2
    P0_PP0 = (P_PP)*(((2+(gamma-1)*M**2)/((gamma+1))))**(gamma/(gamma-1))
    T0_TT0 = (((1+gamma)*M**2*(2+(gamma-1)*M**2))/(1+(gamma*M**2))**2)

    info = messagebox.showinfo('Rayleigh Flow Results', 'For a Mach number of ' + str(M) + ' and a specific heat ratio of ' + str(gamma) + ':\n' + '\n' + 'P/P* = ' + str(P_PP) +'\n'+ 'P0/P0* = ' + str(P0_PP0) + '\n' + 'T/T* = ' + str(T_TT) + '\n' + 'T0/T0* = ' + str(T0_TT0) , parent=ROOT)
    

ROOT = tk.Tk()
ROOT.withdraw()
ans = simpledialog.askstring(title=" Gas Dynamics Calculator",prompt=" Which flow relations are you using? Please enter one of the following: Isentropic, Rayleigh, Normal Shock, or Oblique Shock")

if ans == "Normal Shock":
            M = simpledialog.askfloat("Normal Shock", "Enter the Mach Number : ",parent=ROOT,minvalue=0.0, maxvalue=100000.0)
            gamma = simpledialog.askfloat("Normal Shock", "Enter the specific heat ratio : ",parent=ROOT,minvalue=0.0, maxvalue=100000.0)

            Normal(M, gamma)
                        
elif ans == "normal shock":
            M = simpledialog.askfloat("Normal Shock", "Enter the Mach Number : ",parent=ROOT,minvalue=0.0, maxvalue=100000.0)
            gamma = simpledialog.askfloat("Normal Shock", "Enter the specific heat ratio : ",parent=ROOT,minvalue=0.0, maxvalue=100000.0)

            Normal(M,gamma)
            
elif ans == "Isentropic":
    
            M = simpledialog.askfloat("Isentropic Flow", "Enter the Mach Number : ",parent=ROOT,minvalue=0.0, maxvalue=100000.0)
            gamma = simpledialog.askfloat("Isentropic Flow", "Enter the specific heat ratio : ",parent=ROOT,minvalue=0.0, maxvalue=100000.0)

            Isentropic(M, gamma)

elif ans == "isentropic":
            M = simpledialog.askfloat("Isentropic Flow", "Enter the Mach Number : ",parent=ROOT,minvalue=0.0, maxvalue=100000.0)
            gamma = simpledialog.askfloat("Isentropic Flow", "Enter the specific heat ratio : ",parent=ROOT,minvalue=0.0, maxvalue=100000.0)

            Isentropic(M,gamma)
            
                        
elif ans == "Oblique Shock":
            M = simpledialog.askfloat("Oblique Shock", "Enter the Mach Number : ",parent=ROOT,minvalue=0.0, maxvalue=100000.0)
            gamma = simpledialog.askfloat("Oblique Shock", "Enter the specific heat ratio : ",parent=ROOT,minvalue=0.0, maxvalue=100000.0)
            theta = simpledialog.askfloat("Oblqiue Shock", "Enter the turn angle in degrees : ",parent=ROOT,minvalue=0.0, maxvalue=100000.0)

            Oblique(M,gamma,theta)
            
elif ans == "oblique shock":
            M = simpledialog.askfloat("Oblique Shock", "Enter the Mach Number : ",parent=ROOT,minvalue=0.0, maxvalue=100000.0)
            gamma = simpledialog.askfloat("Oblique Shock", "Enter the specific heat ratio : ",parent=ROOT,minvalue=0.0, maxvalue=100000.0)
            theta = simpledialog.askfloat("Oblqiue Shock", "Enter the turn angle in degrees : ",parent=ROOT,minvalue=0.0, maxvalue=100000.0)

            Oblique(M,gamma,theta)
            
elif ans == "Rayleigh":
    
            M = simpledialog.askfloat("Rayleigh Flow", "Enter the Mach Number : ",parent=ROOT,minvalue=0.0, maxvalue=100000.0)
            gamma = simpledialog.askfloat("Rayleigh Flow", "Enter the specific heat ratio : ",parent=ROOT,minvalue=0.0, maxvalue=100000.0)

            Rayleigh(M, gamma)

elif ans == "rayleigh":
            M = simpledialog.askfloat("Rayleigh Flow", "Enter the Mach Number : ",parent=ROOT,minvalue=0.0, maxvalue=100000.0)
            gamma = simpledialog.askfloat("Rayleigh Flow", "Enter the specific heat ratio : ",parent=ROOT,minvalue=0.0, maxvalue=100000.0)

            Rayleigh(M,gamma)
            
else:
        messagebox.showinfo("Result", "Error, please enter Isentropic, Rayleigh, Normal Shock, or Oblique Shock")


