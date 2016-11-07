!/usr/bin/env python3

# Kuijper, B & Hoyle, R. A.
# The evolution of a maternal effect
# numerical iterations 
#
# This work is licensed under a Creative Commons 
# Attribution-NonCommercial-ShareAlike 4.0 
# International License. 
# http://creativecommons.org/licenses/by-nc-sa/4.0/
#
# when using or extending the code, please cite
# Kuijper, B & Hoyle, R.
# When to rely on maternal effects and when to rely on phenotypic plasticity?
# Evolution, http://dx.doi.org/10.1111/evo.12635 


import os, random, math, copy, sys
import numpy as np
import pandas as pd
from datetime import date, datetime


class IterM:

    # skip in data output
    skip = 1

    # stochastic autocorr
    rho = 0.5

    # mean of stationary autocorrelated noise
    xi_bar = 0; 

    # standard deviation of stationary autocorrelated noise
    sig_xi = 0.01

    # strength of selection on z
    omega_z_2 = 40

    # variance in environmental component
    sigma_e2 = 1.0 

    # intercept of reference envt
    A=0;
    B=2;

    ampl = 1

    # maximal fitness
    Wmax = 1.0

    # initial value of the envt
    epsinit = 0

    # final value of the environment
    epsend = 10

    # length of simulation
#    nrp = 500000
#
#    Uton = 50000
#
#    Utoff = 5000000
    nrp = 50
    Uton = 45
    Utoff = 50

    data_columns = [ 
                    "abar", 
                    "bbar",
                    "mbar",
                    "zbar",
                    "zsbar",
                    "gaz",
                    "gbz",
                    "gmzstar",
                    "fitness",
                    "epst",
                    "epstT",
                    "sigma_z2",
                    "xi_devel",
                    "xi_select"]

    # initialize a single simulation run with particular values
    # for Gaa, Gbb, Gmm, omega_b_2, omega_m_2 which override the class
    # defaults specified above
    # shift: logical variable stating whether we are interested in a single environmental
    # shift (True) or in a sinusoidal environment (False)
    # freq: frequency of the sinusoidal
    # tau: the delay in environmental information
    # 

    def __init__(self, Gaa, Gbb, Gmm, omega_b_2, omega_m_2, shift=True, freq=0, tau=0.25, delta=10, init_g=0.001,init_b=0.001,init_m=0.001):

        self.Gaa = Gaa
        self.Gbb = Gbb
        self.Gmm = Gmm
        self.omega_b_2 = omega_b_2
        self.omega_m_2 = omega_m_2
        self.freq = freq
        self.tau = tau

        # determine magnitude of environmental shift
        # relative to baseline environment
        self.delta = delta

        # initial values of everything
        self.init_m = init_m
        self.init_g = init_g
        self.init_b = init_b

        # print parameters to the screen before running the thing
        print("\n\nGaa:;" + str(self.Gaa) + "\n" + \
                "Gbb:;" + str(self.Gbb) + "\n" + \
                "Gmm:;" + str(self.Gmm) + "\n" + \
                "delta:;" + str(self.delta) + "\n" + \
                "B:;" + str(self.B) + "\n" + \
                "init_g:;" + str(self.init_g) + "\n" + \
                "init_b:;" + str(self.init_b) + "\n" + \
                "init_m:;" + str(self.init_m) + "\n" + \
                "omega_b_2:;" + str(self.omega_b_2) + "\n" + \
                "omega_m_2:;" + str(self.omega_m_2) + "\n" + \
                "omega_z_2:;" + str(self.omega_z_2) + "\n" + \
                "freq:;" + str(self.freq) + "\n" + \
                "tau:;" + str(self.tau))

        # storage for variables
        data = pd.DataFrame(np.nan, 
                index=range(0,10), 
                columns = self.data_columns
        )

        # initialize the first three rows of the relevant
        # columns: mbar, abar, bbar
        data["mbar"][range(0,3)] = self.init_m
        data["abar"][range(0,3)] = self.init_g
        data["bbar"][range(0,3)] = self.init_b
        data["zsbar"][range(0,2)] = 0.0
        data["zbar"][range(0,2)] = 0.0
        data["gmzstar"][range(0,2)] = 0.0
        data["gaz"][range(0,2)] = self.Gaa / (2 - data["mbar"][0])
        data["gbz"][range(0,2)] = self.Gbb * self.epsinit / (2 - data["mbar"][0])
        data["epstT"][range(0,2)] = self.epsinit
        data["epst"][range(0,2)] = self.epsinit

        # initialize phenotypic variance
        data["sigma_z2"][range(0,3)] = (self.Gaa + self.Gbb * self.epsinit**2.0 \
                                                + self.sigma_e2 + 2 * data["mbar"][0] * data["gaz"][0] \
                                                + self.epsinit * data["gbz"][0] + data["gmzstar"]**2 \
                                                + 2 * data["gmzstar"][0] * data["mbar"][0] * data["zsbar"][0] \
                                                + self.Gmm * data["zsbar"]**2)/(1.0 - self.Gmm - data["mbar"][0]**2);

        # create the G-matrix
        self.GG = np.array([[ self.Gaa, 0, 0], [ 0, self.Gbb, 0] , [ 0, 0, self.Gmm]])

        # set the type of environmental change 
        # for the current simulation
        #       -either we have a sudden envt'al shift
        #       -or we have a fluctuating envt
        self.envt_fun = self.envt_shift
        self.envtT_fun = self.envtT_shift

        if not shift:
            self.envt_fun = self.envt_sin
            self.envtT_fun = self.envtT_sin

        # store for the write out
        self.shift = shift

        # assign data to class
        self.data = data
        
        self.init_output()
        self.run()
        self.exit_output()

    # environment changes in a single shift
    # t: the current time t
    def envt_shift(self, t):
        return(self.epsinit + self.delta*(t>self.Uton)*(t<self.Utoff) + (self.epsend - self.epsinit) * (t>=self.Utoff))

    def envtT_shift(self, t):
        return(self.epsinit + self.delta*((t-self.tau)>self.Uton)*((t-self.tau)<self.Utoff) + (self.epsend - self.epsinit) * ((t-self.tau)>=self.Utoff))

    # environment changes periodically 
    # t: the current time t
    def envt_sin(self, t):
        return(self.epsinit+self.ampl*math.sin(self.freq *t))

    def envtT_sin(self, t):
        return(self.epsinit+self.ampl*math.sin(self.freq * (t-self.tau)))

    # run the actual simulation
    def run(self):

        # generate environmental sequence
        self.data["xi_devel"][0] = self.sig_xi * random.gauss(0, 1);

        for k in range(0,self.nrp-1):
            if (self.rho == 0.0):
                self.data["xi_select"][k] = self.sig_xi * random.gauss(0,1);
                self.data["xi_devel"][k] = self.sig_xi * random.gauss(0,1);
            else:
                tmp = (self.rho*self.data["xi_devel"][k]) + (self.sig_xi * math.sqrt(1.0-(self.rho**2))*random.gauss(0,1))
                self.data["xi_select"][k]=tmp

                if (self.tau > 0):
                    j=1.0/self.tau

                    # get through timesteps until you get the actual developmental environment
                    for i in range(0,int(round(j))):
                        tmp =(self.rho*tmp) + (self.sig_xi * math.sqrt(1-(self.rho**2))*random.gauss(0,1))

                    self.data["xi_devel"][k+1] = (self.rho*tmp) + (self.sig_xi * math.sqrt(1-(self.rho**2))*random.gauss(0,1))

        # the guts of the code: loop through time
        for t in range(2,self.nrp-2):
            
            print(t)

            # environmental change
            epsfun = self.envt_fun(t)
            eps_t = epsfun + self.data["xi_select"][t]
            self.data["epst"][t] = eps_t 
            
            epstTfun = self.envtT_fun(t)
            epstT = epstTfun + self.data["xi_devel"][t]
            self.data["epstT"][t] = epstT

            epstT1 = self.data["epstT"][t-1] + self.data["xi_devel"][t-1]
            epstT2 = self.data["epstT"][t-2] + self.data["xi_devel"][t-2]

            # calculate sigma_z^2
            sigma_z2 = self.Gaa + self.Gbb * epstT**2.0 + self.sigma_e2 \
                    + 2 * self.data["mbar"][t] * (self.data["gaz"][t-1] + epstT * self.data["gbz"][t-1])

            sigma_z2 += self.data["gmzstar"][t-1]**2 \
                    + 2*self.data["gmzstar"][t-1]*self.data["mbar"][t]*self.data["zsbar"][t-1] \
                    + self.Gmm * self.data["zsbar"][t-1]**2

            sigma_z2 += (self.Gmm + self.data["mbar"][t]**2) * self.data["sigma_z2"][t-1];

            sigma_z2 = max(sigma_z2,0);

#            sigma_z2 = sigma_z2 / (1.0 - self.Gmm - self.data["mbar"][t]**2)

            self.data["sigma_z2"][t] = sigma_z2

            # calculate gamma
            gamma = 1.0 / (self.omega_z_2 + sigma_z2)
            gamma_b = 1.0 / (self.omega_b_2 + self.Gbb)
            gamma_m = 1.0 / (self.omega_m_2 + self.Gmm)

            # calculate changes in elevation and plasticity
            # update \bar{z}_{t}
            self.data["zbar"][t] = self.data["abar"][t] + (self.data["bbar"][t]*epstT) + self.data["gmzstar"][t-1] + self.data["mbar"][t]*self.data["zsbar"][t-1]

            # update the value of theta
            theta = self.A + self.B * eps_t

            tmp = self.data["zbar"][t] - theta
            mat1 = np.multiply(-1.0/self.omega_z_2,np.array([
                    (1+self.data["mbar"][t])*tmp,
                    (epstT+self.data["mbar"][t]*epstT1)*tmp,
                    (self.data["zsbar"][t-1]+self.data["mbar"][t]*self.data["zsbar"][t-2]+0.5*self.data["gmzstar"][t-1])*tmp
                    ]))

            dsiga = 2 * (self.data["mbar"][t] * self.data["gmzstar"][t-1] + self.Gmm * self.data["zsbar"][t-1])/(1.0 - self.Gmm - self.data["mbar"][t]**2)

            dsigb = dsiga * epstT1

            dsigm = (2.0 / (1.0 - self.Gmm - self.data["mbar"][t]**2)) * (self.data["gaz"][t-1] * (1.0 + 0.5 * self.data["mbar"][t]) + (epstT + 0.5 * self.data["mbar"][t] * epstT1) * self.data["gbz"][t-1] + 0.5 * self.data["gmzstar"][t-1]**2 + self.data["gmzstar"][t-1] * self.data["zsbar"][t-1] * (1.0 + 0.5 * self.data["mbar"][t]) + self.data["gmzstar"][t-1] * self.data["mbar"][t] * self.data["zsbar"][t-2] + self.Gmm * self.data["zsbar"][t-1] * self.data["zsbar"][t-2]) + 2 * self.data["mbar"][t] * sigma_z2 / (1-self.Gmm-self.data["mbar"][t]**2)

            mat2 = np.multiply(0.5 * (-1.0 / self.omega_z_2),np.array([ dsiga, dsigb, dsigm ]))

            mat3 = np.multiply(- (1.0 / self.omega_z_2), np.array([ 0, self.omega_z_2 * self.data["bbar"][t] / self.omega_b_2, self.omega_z_2 * self.data["mbar"][t] / self.omega_m_2 ]))

            beta = mat1 + mat2 + mat3

            delta_change = np.dot(self.GG,beta)

            prefacfit = math.sqrt(gamma * gamma_b * gamma_m * self.omega_z_2 * self.omega_b_2 * self.omega_m_2)

            self.data["fitness"][t] = self.Wmax * prefacfit * math.exp(-(gamma/2.0) * ((self.data["zbar"][t] - theta)**2) - .5 * gamma_m * (self.data["mbar"][t]**2) - .5 * gamma_b * (self.data["bbar"][t]**2))

            self.data["abar"][t+1] = self.data["abar"][t] + delta_change[0]
            self.data["bbar"][t+1] = self.data["bbar"][t] + delta_change[1]
            self.data["mbar"][t+1] = self.data["mbar"][t] + delta_change[2]

            self.data["zsbar"][t] = self.data["abar"][t+1] + self.data["bbar"][t+1] * epstT + self.data["gmzstar"][t-1] + self.data["mbar"][t+1] * self.data["zsbar"][t-1]

            # update C_{m_{t}z_{t-1}^*}
            self.data["gmzstar"][t] = 0.5 * self.data["mbar"][t+1] * self.data["gmzstar"][t-1] + 0.5 * self.Gmm * self.data["zsbar"][t-1]

            self.data["gaz"][t] = 0.5 * self.Gaa + 0.5 * self.data["mbar"][t+1] * self.data["gaz"][t-1]

            self.data["gbz"][t] = 0.5 * self.Gbb * epstT + 0.5 * self.data["mbar"][t+1] * self.data["gbz"][t-1]

            # get last index
            self.update_dataframe(t)

    def update_dataframe(self,t):
            
        idxmax = self.data["abar"].index[-1]
        idxmin = self.data["abar"].index[0]
        self.data = self.data.append(pd.DataFrame(np.nan, index=[idxmax+1], columns=self.data_columns))

        if t % self.skip == 0:
            self.output_file.write(";".join([str(x) for x in self.data.ix[t]])+"\n")

        self.data = self.data.drop(idxmin)
        print(self.data)


    def init_output(self):

        d = datetime.now()
        filename = "iter_" + d.strftime("%d_%m_%Y_%H%M%S_%f") + ".csv"

        self.output_file = open(filename,"w")

        self.output_file.write(";".join(self.data_columns) + "\n")

    def exit_output(self):
        self.output_file.write("\n\ntype:;cascading_m\n" + \
                "Gaa:;" + str(self.Gaa) + "\n" + \
                "Gbb:;" + str(self.Gbb) + "\n" + \
                "Gmm:;" + str(self.Gmm) + "\n" + \
                "delta:;" + str(self.delta) + "\n" + \
                "B:;" + str(self.B) + "\n" + \
                "A:;" + str(self.A) + "\n" + \
                "omega_b_2:;" + str(self.omega_b_2) + "\n" + \
                "omega_m_2:;" + str(self.omega_m_2) + "\n" + \
                "omega_z_2:;" + str(self.omega_z_2) + "\n" + \
                "shift:;" + str(self.shift) + "\n" + \
                "sig_xi:;" + str(self.sig_xi) + "\n" + \
                "rho:;" + str(self.rho) + "\n" + \
                "freq:;" + str(self.freq) + "\n" + \
                "tau:;" + str(self.tau)) 

        self.output_file.close()

# end class 

Gaa = float(sys.argv[1])
Gbb = float(sys.argv[2])
Gmm = float(sys.argv[3])
omega_b_2 = float(sys.argv[4])
omega_m_2 = float(sys.argv[5])
shift = bool(int(sys.argv[6]))
freq = float(sys.argv[7])
tau = float(sys.argv[8])
delta = float(sys.argv[9])
init_g = float(sys.argv[10])
init_b = float(sys.argv[11])
init_m = float(sys.argv[12])

a = IterM(Gaa, Gbb, Gmm, omega_b_2, omega_m_2, shift, freq, tau, delta, init_g, init_b, init_m)
