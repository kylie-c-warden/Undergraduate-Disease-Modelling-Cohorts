import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from SALib.sample import saltelli
from SALib.analyze import sobol

# Define the ODE system
def odes(x, t, u_c, u_p, u_i, S_a, S_c, S_p, S_i, alpha_c, Beta, r, lambda_p, lambda_i):
    A, P, C, I = x  # State variables

    # ODE system
    dPdt = lambda_p - S_p * P - Beta * C * P
    dCdt = u_c * (1 - alpha_c) * C - S_c * (1 + alpha_c) * I * C + Beta * C * P
    dAdt = S_p * P + S_c * (1 + alpha_c) * I * C - r * I * A
    dIdt = lambda_i - S_i * I + u_c * C

    return [dAdt, dPdt, dCdt, dIdt]

# Define the problem and parameter ranges
problem = {
    "num_vars": 12,
    "names": ["u_c", "u_p", "u_i", "S_a", "S_c", "S_p", "S_i", "alpha_c", "Beta", "r", "lambda_p", "lambda_i"],
    "bounds": [
        [0.02, 0.05],  # u_c
        [0.01, 0.05],  # u_p
        [0.01, 0.05],  # u_i
        [0.01, 0.04], # S_a
        [0.009, 0.03], # S_c
        [0.01, 0.05], # S_p
        [0.005, 0.02], # S_i
        [0, 0.2],      # alpha_c
        [0.01, 0.3],   # Beta
        [0.01, 0.05],  # r
        [0.01, 0.05],  # lambda_p
        [0.01, 0.8]    # lambda_i
    ]
}

# Generate Sobol samples
param_values = saltelli.sample(problem, 1024)  # Increase samples for higher accuracy

# Initial conditions and time range
x_0 = [4000, 1.75e6, 40000, 10e6]
t = np.linspace(0, 100, 1000)  # Simulation time

# Run simulations for each sample
outputs_A = np.zeros((param_values.shape[0], len(t)))
outputs_P = np.zeros((param_values.shape[0], len(t)))
outputs_C = np.zeros((param_values.shape[0], len(t)))
outputs_I = np.zeros((param_values.shape[0], len(t)))

for i, params in enumerate(param_values):
    sol = odeint(odes, x_0, t, args=tuple(params))  # Solve ODEs
    outputs_A[i] = sol[:, 0]  # Extract A(t)
    outputs_P[i] = sol[:, 1]  # Extract P(t)
    outputs_C[i] = sol[:, 2]  # Extract C(t)
    outputs_I[i] = sol[:, 3]  # Extract I(t)

# Compute Sobol indices for each variable at the final time point (t = 100)
Y_A = outputs_A[:, -1]
Y_P = outputs_P[:, -1]
Y_C = outputs_C[:, -1]
Y_I = outputs_I[:, -1]

Si_A = sobol.analyze(problem, Y_A)
Si_P = sobol.analyze(problem, Y_P)
Si_C = sobol.analyze(problem, Y_C)
Si_I = sobol.analyze(problem, Y_I)

# Function to plot Sobol indices
def plot_sobol_indices(Si, title):
    plt.figure(figsize=(10, 6))
    plt.bar(problem["names"], Si["ST"], color="blue", alpha=0.7, label="Total-Effect Index")
    plt.bar(problem["names"], Si["S1"], color="red", alpha=0.7, label="First-Order Index")
    plt.xlabel("Parameters")
    plt.ylabel("Sobol Sensitivity Index")
    plt.title(title)
    plt.legend()
    plt.xticks(rotation=45)
    plt.grid()
    plt.show()

# Plot Sobol indices for each output
plot_sobol_indices(Si_A, "Sobol Sensitivity Analysis for A(t) at t=100")
plot_sobol_indices(Si_P, "Sobol Sensitivity Analysis for P(t) at t=100")
plot_sobol_indices(Si_C, "Sobol Sensitivity Analysis for C(t) at t=100")
plot_sobol_indices(Si_I, "Sobol Sensitivity Analysis for I(t) at t=100")

# Print Sobol indices for reference
print("Sobol Indices for A(t) at t=100")
print("First-order (S1):", Si_A["S1"])
print("Total-effect (ST):", Si_A["ST"])

print("\nSobol Indices for P(t) at t=100")
print("First-order (S1):", Si_P["S1"])
print("Total-effect (ST):", Si_P["ST"])

print("\nSobol Indices for C(t) at t=100")
print("First-order (S1):", Si_C["S1"])
print("Total-effect (ST):", Si_C["ST"])

print("\nSobol Indices for I(t) at t=100")
print("First-order (S1):", Si_I["S1"])
print("Total-effect (ST):", Si_I["ST"])
