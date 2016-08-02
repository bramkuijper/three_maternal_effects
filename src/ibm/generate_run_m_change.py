#!/usr/bin/env python

import numpy, math

mu_g = 0.01

mu_m_combis = [
                [ 0.01, 0, 0],
                [ 0, 0.01, 0],
                [ 0, 0, 0.01]
                ]

mu_b = [0]

sdmu = 0.02

sigma_e = 1.0

sigma_ksi = 0

wmin = 0

rho_t = 0.5

omega2 = [ 0.7 ]

tau = [ 0.25 ]

rate1 = list(numpy.arange(0,math.pi+math.pi/10,math.pi/10))
int_ptb = list(numpy.arange(0,5,1.0))
ampl1 = 1.0

exe = "./xm_sin_change"

reps = 5

ctr = 0

for rep_i in range(0,reps):
    for mu_m_combi_i in mu_m_combis:
        for mu_b_i in mu_b:
            for omega2_i in omega2: 
                for tau_i in tau:
                    for rate1_i in rate1: 
                        for int_ptb_i in int_ptb:
                            
                            print("echo " + str(ctr))

                            ctr+=1
                            
                            print(exe + " " + str(mu_g) + "\t" + str(mu_m_combi_i[0]) + " " + str(mu_m_combi_i[1]) + " " + str(mu_m_combi_i[2]) + "\t" + str(mu_b_i) + \
                                    " " + str(sdmu) + "\t" + str(sigma_e) + " " + str(sigma_ksi) + "\t" + str(wmin) + " " + str(omega2_i) + \
                                    " " + " 100 100 100 100 100 " + str(tau_i) + " " + " 0 " + str(rate1_i) + " 1.0 " + str(int_ptb_i) + " " + str(rate1_i) + " 1.0")


