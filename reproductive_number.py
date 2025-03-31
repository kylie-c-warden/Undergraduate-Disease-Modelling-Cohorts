alpha_c = 0    #keep curcumin at a range of 0.1 to 0.2

Beta = 0.0000012     # 0.00012   # 0.0001597

# Deltas
S_c = 0.9      # 0.9    # 0.015 
S_p = 9e-6     # 9e-6   # 0.00399995 
S_i = 2    #3.162   # 0.0527

lambda_p = 18       #18       # 0.3
lambda_i = 0.3      #0.3      # 0.005

# mu 
u_c = 1.8      # 1.8    # 0.03 per hr 

r = 1.8             # 1.8     # 0.03

#reproductive number
r_0 = (u_c * (1 - alpha_c) + Beta * (lambda_p / S_p)) / ((S_c * (1 + alpha_c)) * (lambda_i / S_i)) 


