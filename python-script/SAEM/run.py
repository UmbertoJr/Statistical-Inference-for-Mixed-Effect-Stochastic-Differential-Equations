from SAEM import Stochasti_Approximation_Expectation_Maximization
from time import time

start = time()
# self.PARMIN = [0.01, 1e-6, 0.0001, 1e-5, 0.0001, 1e-5] ; self.PARMAX = [1, 1e-1, 0.05, 1e-1, 1, 1e-3]
theta = [0.3, 0.009, 0.2, 0.3, 0.3, 0.0001]   # [0.27, 0.0005, 0.04, 0.0035, 0.18, 0.0001]
SAEM = Stochasti_Approximation_Expectation_Maximization()
#SAEM.inizializza(theta)
SAEM.run(theta)
print("tempo esecuzione : %.5f "% (time()-start))

import matplotlib.pyplot as plt
import numpy as np
def x_fun_espl(t, alpha, kappa, Y_0):
    y = kappa * np.exp( np.log(Y_0/kappa) *  np.exp(- alpha* t))
    return y

s_= 1
SAEM.select_subj(s_)
plt.plot(SAEM.timei, SAEM.yobsi, "r*")
tempo = np.linspace(0.0, 25, 25000)
x_hat_medio = x_fun_espl(tempo,  SAEM.theta_choose[0], SAEM.theta_choose[1], 1e-7)  
x_hati = x_fun_espl(tempo, SAEM.PhiALL[SAEM.iterazione,s_ - 1, 0], SAEM.PhiALL[SAEM.iterazione,s_ -1, 1], 1e-7)  
plt.plot(tempo, x_hati, "g-", label = "ultimi phii") 
plt.plot(tempo, x_hat_medio, label = "medio popolazione")
plt.legend()