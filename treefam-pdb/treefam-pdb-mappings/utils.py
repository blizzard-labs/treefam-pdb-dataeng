import numpy as np
import pandas as pd
import json
from json import JSONEncoder

def numpy2csv(array, name):
    DF = pd.DataFrame(array)
    DF.to_csv(name)

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

def cleanUp(keySegs_array, minL):
    newArray = []
    for segment in keySegs_array:
        if len(segment) != 0 and len(segment) > minL:
            newArray.append(str(segment[0] + 1) + '-' + str(segment[-1] + 1))
    return newArray
    