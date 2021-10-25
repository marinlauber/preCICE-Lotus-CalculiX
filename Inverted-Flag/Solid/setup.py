import numpy as np
import matplotlib.pyplot as plt

def get_rho_s(Mratio, L, h):
    return Mratio*L/h

def get_E(kappa, L, h, nu=0.4):
    return kappa*L**3*12*(1-nu**2)/(h**3)

L = 64
h = 2.56

print(get_rho_s(0.5, L, h))
print(get_E(0.4, L, h))
