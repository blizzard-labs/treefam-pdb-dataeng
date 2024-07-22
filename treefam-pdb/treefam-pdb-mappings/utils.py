import numpy as np
import pandas as pd
import json
from json import JSONEncoder

#Convert NumPy to CSV
def numpy2csv(array, name):
    DF = pd.DataFrame(array)
    DF.to_csv(name)

#Convert CSV to NumPy
def csv2numpy(csv):
    raw = np.genfromtxt(csv, delimiter=',')
    return (raw[1:raw.shape[0], 1:raw.shape[1]])

class NumpyArrayEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return JSONEncoder.default(self, obj)

def numpy2json(npArray, jFile):
    data = {'array': npArray}
    with open(jFile, "w") as write_file:
        json.dump(data, write_file, cls=NumpyArrayEncoder)

def json2numpy(jFile):
    with open(jFile, "r") as read_file:
        decoded = json.load(read_file)
        return(np.asarray(decoded['array']))

if __name__ == '__main__':
    numpy2json(np.array([[11 ,22, 33], [44, 55, 66], [77, 88, 99]]), "mapping.json")
    print(json2numpy("mapping.json"))
    
    