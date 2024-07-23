import numpy as np

anarray = np.zeros((4, 4))
anotherarray = np.zeros((3, 1))

alist = [anarray, anotherarray]
print(alist[0].shape)