import pandas as pd
import numpy as np
import random


ax = pd.DataFrame(np.nan,index=range(0,10),columns=["N","K"])

ax["N"][0] = 100

K = 500

r = 1.01


for t in range(1,50):

    ntmin1 = ax["N"][t-1]

    ax["N"][t] = r * ntmin1 * (1.0 - ntmin1/K)

    ax["K"][t] = random.gauss(0,1)

    if t > 8:
        idxmax = ax["N"].index[-1]
        ax = ax.append(pd.DataFrame(np.nan,index=[idxmax+1],columns=["N","K"]))
        idxmin = ax["N"].index[0]

        print(ax)
        print(t)
        print(";".join([str(x) for x in ax.ix[t]]))

        ax = ax.drop(idxmin)



k
