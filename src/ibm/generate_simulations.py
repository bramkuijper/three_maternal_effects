#!/usr/bin/env python3

import numpy as np
import math

mu_g = 0.01
wmin = 0
omega2 = [ 0.7 ]
mu_m_g = [ 0, 0.01 ]
mu_m_e = [ 0, 0.01 ]
mu_m_m = [ 0, 0.01 ]
mu_b = [ 0, 0.01 ]
sigma_e = [ 1.0 ]
sigma_ksi = 0.01

step = 30
freq = list(np.arange(0, math.pi + math.pi/step, math.pi/step))

tau = 0.5
sdmu = 0.02

replicates = 10
ctr = 0 
exe = "./xm_sin"

for rep_i in range(0, replicates):
    for mu_m_g_i in mu_m_g:
        for mu_m_e_i in mu_m_e:
            for mu_m_m_i in mu_m_m:
                for mu_b_i in mu_b:
                    for sigma_e_i in sigma_e:
                        for freq_i in freq:
                            for omega2_i in omega2:

                                print("echo " + str(ctr))
                                ctr+=1

                                print(exe + " " 
                                        + str(mu_g) + " " 
                                        + str(mu_m_g_i) + " "
                                        + str(mu_m_e_i) + " "
                                        + str(mu_m_m_i) + " "
                                        + str(mu_b_i) + " "
                                        + str(sdmu) + " "
                                        + str(freq_i) + " "
                                        + str(sigma_e_i) + " "
                                        + str(sigma_ksi) + " "
                                        + str(wmin) + " "
                                        + str(0) + " "
                                        + str(omega2_i) + " "
                                        + str(100) + " "
                                        + str(100) + " "
                                        + str(100) + " "
                                        + str(100) + " "
                                        + str(tau) + " ")

