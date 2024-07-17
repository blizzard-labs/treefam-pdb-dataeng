import numpy as np

new = np.zeros((4, 4))
print(new)

def optional(something, **abc):
    if abc:
        print(abc)
        print(np.array(abc["abc"]).shape)
    else: print("hi")

optional(1)   