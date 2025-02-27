#import matplotlib.pyplot as plt
#from scipy.integrate import odeint
import numpy as np

#mu =
u_c = 0.03
u_p = 0.03
u_i = 0.033

#delta = S
S_a = 0.02 
S_c = 0.015   
S_p = 0.03  
S_i = 0.01

alpha_c = 0.2   
Beta = 0.20  
r = 0.03
lambda_p = 0.03
lambda_i = 0.05

r_0 = (u_c * (1 - alpha_c) + Beta * (lambda_p / S_p)) / ((S_c * (1 + alpha_c)) * (lambda_i / S_i)) 

#print(r_0)

