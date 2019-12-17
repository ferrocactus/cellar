import numpy as np

def read_data():
    return np.random.random(size=(100, 50))

if __name__ == '__main__':
    print(read_data())