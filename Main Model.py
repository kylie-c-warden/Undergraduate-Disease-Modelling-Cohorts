
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import odeint
import numpy as np
from reproductive_number import *



def odes (x, t): 

    A = x[0]
    P = x[1]
    C = x[2]
    I = x[3]

    # our ODEs 
    dPdt = lambda_p - S_p * P - Beta * C * P   #good 
    dCdt = u_c * (1 - alpha_c) * C  - S_c * (1+ alpha_c) * I * C + Beta * C * P #good
    dAdt = S_p * P + S_c * (1 + alpha_c) * I * C - r * I * A #good
    dIdt = lambda_i - S_i * I + u_c * C #good 

    return [dAdt, dPdt, dCdt, dIdt]

#Sensitivity of parameters:
S_alphaC = (alpha_c / r_0)* ((S_i*(-1)*(2*S_p*u_c - Beta*lambda_p)) / (S_p*S_c*lambda_i*(1+alpha_c)**2))

S_beta = (Beta / r_0)*((lambda_p*S_i)/(S_p*S_c*(1+alpha_c)*lambda_i))

S_deltaC = (S_c / r_0)*(-1)*(S_i*(u_c*(1-alpha_c)*S_p+Beta*lambda_p)/(((S_c)**2)*(1+alpha_c)*lambda_i*S_p))

S_deltaI = (S_i / r_0)*((u_c*(1-alpha_c)*S_p+Beta*lambda_p)/(S_c*S_p*(1+alpha_c)*lambda_i))

S_deltaP = (S_p / r_0)*(-1)*((Beta*lambda_p*S_i) / ((S_c*(S_p)**2) *(1+alpha_c)*lambda_i))

S_lambdaI = (lambda_i / r_0)*(-1)*S_i*((u_c*(1-alpha_c)+Beta*(lambda_p / S_p))/(S_c*(1+alpha_c)*((lambda_i)**2)))

S_lambdaP = (lambda_p / r_0)*(Beta*S_i/S_p*S_c*(1+alpha_c)*lambda_i)

S_muC = (u_c/r_0)*((S_i*(1-alpha_c))/(S_c*(1+alpha_c)*lambda_i))

params = dict(aC = S_alphaC, B = S_beta, dC = S_deltaC, dI = S_deltaI, dP = S_deltaP, li = S_lambdaI, lP = S_lambdaP, mC = S_muC)
print(params)
#  initial conditions 
# x_0 = [a_0, p_0, c_0, i_0] 

x_0 = [4000, 1.75e6, 40000, 10e6] 


print(odes(x_0,0))
tf = 100
t = np.linspace(0,tf,1000)
x = odeint(odes, x_0, t)

A = x[:,0]
P = x[:,1]
C = x[:,2]
I = x[:,3]


# # Plotting each equation in a separate subplot with a log y-axis
fig, axis = plt.subplots(2, 3, figsize=(14, 7))

# Plot P(t)
axis[0, 0].plot(t, P, color='blue')
axis[0, 0].set_title("$P(t)$ Graph")
# axis[0, 0].set_yscale('log')
axis[0, 0].set_xlabel('$t$')
axis[0, 0].set_ylabel('$P(t)$')
axis[0, 0].grid()

# Plot A(t)
axis[0, 1].plot(t, A, color='green')
axis[0, 1].set_title("$A(t)$ Graph")
# axis[0, 1].set_yscale('log')
axis[0, 1].set_xlabel('$t$')
axis[0, 1].set_ylabel('$A(t)$')
axis[0, 1].grid()

# Plot C(t)
axis[1, 0].plot(t, C, color='red')
axis[1, 0].set_title("$C(t)$ Graph")
# axis[1, 0].set_yscale('log')
axis[1, 0].set_xlabel('$t$')
axis[1, 0].set_ylabel('$C(t)$')
axis[1, 0].grid()

# Turn off the last empty subplot
axis[1,1].plot(t, I, color='purple')
axis[1,1].set_title('$I(t)$ Graph')
axis[1, 1].set_xlabel('$t$')
axis[1,1].set_ylabel('$I(t)$')
axis[1,1].grid()

#sensitivity
x = params.keys()
y = params.values()
axis[0,2].bar(x, y, width=.5)

#make last spot empty
axis[1,2].axis("off")

# Adjust layout for better spacing
plt.tight_layout()
plt.show()