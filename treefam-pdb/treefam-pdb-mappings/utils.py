import numpy as np
import pandas as pd

#Convert NumPy to CSV
def numpy2csv(array, name):
    DF = pd.DataFrame(array)
    DF.to_csv(name)

#Convert CSV to NumPy
def csv2numpy(csv):
    raw = np.genfromtxt(csv, delimiter=',')
    return (raw[1:raw.shape[0], 1:raw.shape[1]])