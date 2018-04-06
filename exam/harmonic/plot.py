import os
join = os.path.join
exists = os.path.exists
expanduser = os.path.expanduser
import sys
import matplotlib as mpl
mpl.use('PDF')
import matplotlib.pylab as plt
import numpy as np
tr = np.transpose
import pandas as pd

import mangan4 as mn
import json


def run():
    man = mn.Man()
    x  = man['_x', 0]
    ts = man['t',:]    
    for it in [0, 4, 8]:
        fr = man['fr',it]
        f  = fr[0::2] + 1j*fr[1::2]
        plt.plot(x, f.real, label="calc,t={0}".format(ts[it]))
    plt.legend()
    plt.xlim(-5.0, 5.0)
    plt.savefig("rho.pdf")

if __name__=="__main__":
    run()

    
