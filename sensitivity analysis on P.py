import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define the system of ODEs
def odes(x, t, params):
    # Unpack parameters
    u_c, u_p, u_i, S_a, S_c, S_p, S_i, alpha_c, Beta, r, lambda_p, lambda_i = params
    
    # State variables
    A, P, C, I = x

    # ODE system
    dPdt = lambda_p - S_p * P - Beta * C * P
    dCdt = u_c * (1 - alpha_c) * C  - S_c * (1 + alpha_c) * I * C + Beta * C * P
    dAdt = S_p * P + S_c * (1 + alpha_c) * I * C - r * I * A
    dIdt = lambda_i - S_i * I + u_c * C

    return [dAdt, dPdt, dCdt, dIdt]

# Baseline parameters
baseline_params = [0.03, 0.03, 0.033, 0.02, 0.015, 0.03, 0.01, 0, 0.20, 0.03067, 0.03, 0.05]

# Initial conditions
x_0 = [4000, 1.75e6, 40000, 10e6]
tf = 100
t = np.linspace(0, tf, 1000)

# Function to run simulation
def run_simulation(params):
    return odeint(odes, x_0, t, args=(params,))

# Run baseline simulation
baseline_solution = run_simulation(baseline_params)

# Sensitivity Analysis: Vary Each Parameter by ±10%
sensitivity_results = {}
perturbation = 0.1  # 10% variation

for i, param_name in enumerate(["u_c", "u_p", "u_i", "S_a", "S_c", "S_p", "S_i", "alpha_c", "Beta", "r", "lambda_p", "lambda_i"]):
    # Increase parameter by 10%
    params_high = baseline_params.copy()
    params_high[i] *= (1 + perturbation)
    sol_high = run_simulation(params_high)

    # Decrease parameter by 10%
    params_low = baseline_params.copy()
    params_low[i] *= (1 - perturbation)
    sol_low = run_simulation(params_low)

    # Store results
    sensitivity_results[param_name] = {"high": sol_high, "low": sol_low}

# Plot sensitivity results for P(t)
plt.figure(figsize=(10, 6))
plt.plot(t, baseline_solution[:, 1], label="Baseline", color="black")

for param_name, result in sensitivity_results.items():
    plt.plot(t, result["high"][:, 1], '--', label=f"{param_name} +10%")
    plt.plot(t, result["low"][:, 1], '--', label=f"{param_name} -10%")

plt.xlabel("Time (t)")
plt.ylabel("$P(t)$")
plt.title("Sensitivity Analysis on $P(t)$")
plt.legend()
plt.grid()
plt.show()
